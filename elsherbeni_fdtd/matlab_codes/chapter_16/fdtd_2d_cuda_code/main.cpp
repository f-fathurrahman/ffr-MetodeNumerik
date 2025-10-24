#include <stdio.h>
#include <math.h>
#include <time.h>

#include <GL/glew.h>
#include <GL/glu.h>
#include <GL/freeglut.h>

#include <cuda.h>
#include <cudaGL.h>
#include <cuda_runtime_api.h>

#include "global.h"

// global variables
unsigned int window_width = 512;
unsigned int window_height = 512;
unsigned int image_width = 512;
unsigned int image_height = 512;
int iGLUTWindowHandle = 0;          // handle to the GLUT window
size_t number_of_bytes; 

GLuint pbo_destination;
struct cudaGraphicsResource *cuda_pbo_destination_resource;
GLuint cuda_result_texture;	

time_t start, stop;

////////////////////////////////////////////////////////////////////////////////
extern "C" void
createImageOnGpu(unsigned int* g_odata);

// Forward declarations
void Cleanup(int iExitCode);

// GL functionality
bool initializeGL(int argc, char** argv);

void createPixelBufferObject(GLuint* pbo, struct cudaGraphicsResource **pbo_resource);
void deletePBO(GLuint* pbo);

void createTextureDestination(GLuint* cuda_result_texture, unsigned int size_x, unsigned int size_y);
void deleteTexture(GLuint* tex);

// rendering callbacks
void idle();
void keyboard(unsigned char key, int x, int y);
void reshape(int w, int h);

void setImageAndWindowSize()
{
	image_width  = nxx;
	image_height = nyy;

	if (nxx>nyy)
		window_height = window_width*nyy/nxx;
	else
		window_width = window_height*nxx/nyy;
}

////////////////////////////////////////////////////////////////////////////////
void createPixelBufferObject(GLuint* pbo, struct cudaGraphicsResource **pbo_resource)
{
    // set up vertex data parameter
	unsigned int texture_size;

    texture_size = sizeof(GLubyte) * image_width * image_height * 4;
    void *data = malloc(texture_size);

    // create buffer object
    glGenBuffers(1, pbo);
    glBindBuffer(GL_ARRAY_BUFFER, *pbo);
    glBufferData(GL_ARRAY_BUFFER, texture_size, data, GL_DYNAMIC_DRAW);
    free(data);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // register this buffer object with CUDA
    cudaGraphicsGLRegisterBuffer(pbo_resource, *pbo, cudaGraphicsMapFlagsNone);
}

void deletePBO(GLuint* pbo)
{
    glDeleteBuffers(1, pbo);
    *pbo = 0;
}

// runIterationAndDisplay image to the screen as textured quad
void displayTextureImage(GLuint texture)
{
    glBindTexture(GL_TEXTURE_2D, texture);
    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(domain_min_x, domain_min_x+nxx*dx, domain_min_y, domain_min_y+nyy*dy, -1.0, 1.0);

    glMatrixMode( GL_MODELVIEW);
    glViewport(0, 0, window_width, window_height);

    glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0); glVertex3f(domain_min_x, domain_min_y, 0.0);
    glTexCoord2f(1.0, 0.0); glVertex3f(domain_min_x+nxx*dx, domain_min_y, 0.0);
    glTexCoord2f(1.0, 1.0); glVertex3f(domain_min_x+nxx*dx, domain_min_y+nyy*dy, 0.0);
    glTexCoord2f(0.0, 1.0); glVertex3f(domain_min_x, domain_min_y+nyy*dy, 0.0);
    glEnd();

    glDisable(GL_TEXTURE_2D);

	glCallList(objects_display_list);

}

void runIterationAndDisplay()
{
	// run a number of FDTD iterations on GPU using CUDA
	for (int i=0; i< plotting_step; i++)
	{
		if (time_step<number_of_time_steps) 
			fdtdIterationOnGpu();
		else
		{
			time(&stop);
			printf("Time elapsed in loop: %f sec\n ",(double )difftime(stop, start));

			fetchResultsFromGpuMemory();
			deallocateCudaArrays();
			deallocateArrays();
			saveSampledFieldsToFile();
			Cleanup(1);
		}
	}

	// Create image of field using CUDA
	unsigned int* image_data;
	// map the GL buffer to CUDA
	cudaGraphicsMapResources(1, &cuda_pbo_destination_resource, 0);    
	cudaGraphicsResourceGetMappedPointer((void **)&image_data, &number_of_bytes, cuda_pbo_destination_resource);

	// execute CUDA kernel
	createImageOnGpu(image_data);
	// unmap the GL buffer
	cudaGraphicsUnmapResources(1, &cuda_pbo_destination_resource, 0);

	// Create a texture from the buffer
	glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, pbo_destination);
	glBindTexture(GL_TEXTURE_2D, cuda_result_texture);
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, image_width, image_height, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glBindBuffer(GL_PIXEL_PACK_BUFFER_ARB, 0);
	glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);

	// draw the image
 	displayTextureImage(cuda_result_texture);

	cudaThreadSynchronize();

	// swap the front and back buffers
	glutSwapBuffers();
}

void idle()
{
    glutPostRedisplay();
}

////////////////////////////////////////////////////////////////////////////////
void keyboard(unsigned char key, int /*x*/, int /*y*/)
{
    switch(key) {
    case(27) :
        Cleanup(1);
    case ' ':
        break;
    case 'a':
        break;
    case '=':
    case '+':
        break;
    case '-':
        break;
    }
}

void reshape(int w, int h)
{
    window_width = w;
    window_height = h;
}

////////////////////////////////////////////////////////////////////////////////
void createTextureDestination(GLuint* cuda_result_texture, unsigned int size_x, unsigned int size_y)
{
    // create a texture
    glGenTextures(1, cuda_result_texture);
    glBindTexture(GL_TEXTURE_2D, *cuda_result_texture);

    // set basic parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, size_x, size_y, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
}


////////////////////////////////////////////////////////////////////////////////
void deleteTexture(GLuint* tex)
{
    glDeleteTextures(1, tex);
    *tex = 0;
}

////////////////////////////////////////////////////////////////////////////////
void initializeGLBuffers()
{
    // create pixel buffer object
    createPixelBufferObject(&pbo_destination, &cuda_pbo_destination_resource);
    // create texture that will receive the result of CUDA
    createTextureDestination(&cuda_result_texture, image_width, image_height);     
}



////////////////////////////////////////////////////////////////////////////////
void Cleanup(int iExitCode)
{
    cudaGraphicsUnregisterResource(cuda_pbo_destination_resource);
	deletePBO(&pbo_destination);
    deleteTexture(&cuda_result_texture);
    cudaThreadExit();
    if(iGLUTWindowHandle)glutDestroyWindow(iGLUTWindowHandle);
    exit (iExitCode);
}

// Initialize GL
bool initializeGL(int argc, char **argv )
{
	setImageAndWindowSize();

	// Create GL context
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(window_width, window_height);
    iGLUTWindowHandle = glutCreateWindow("CUDA OpenGL FDTD");

    // initialize necessary OpenGL extensions
    glewInit();
    if (! glewIsSupported(
        "GL_VERSION_2_0 " 
        "GL_ARB_pixel_buffer_object "
        "GL_EXT_framebuffer_object "
        )) {
        printf("ERROR: Support for necessary OpenGL extensions missing.");
        fflush(stderr);
        return false;
    }
	
	// Initialize GLUT event functions
	glutDisplayFunc(runIterationAndDisplay);
    glutKeyboardFunc(keyboard);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);

    return true;
}
void createDisplayListForObjects()
{
	float min_x, min_y, max_x, max_y;
	float pi = 4*atan(1.0);
	float phi, x, y;

	objects_display_list = glGenLists(1);
	glNewList(objects_display_list, GL_COMPILE);

	// draw sources
		glColor3f(1.0f, 0.0f, 0.0f);
		for (int i=0;i<number_of_impressed_J;i++)
		{
			min_x = impressed_J_min_x[i];
			min_y = impressed_J_min_y[i];
			max_x = impressed_J_max_x[i];
			max_y = impressed_J_max_y[i];
			glBegin(GL_LINES);
			glVertex3f(min_x, min_y, -0.1f); 
			glVertex3f(min_x, max_y, -0.1f); 
			glVertex3f(min_x, min_y, -0.1f); 
			glVertex3f(max_x, min_y, -0.1f); 

			glVertex3f(max_x, min_y, -0.1f); 
			glVertex3f(max_x, max_y, -0.1f); 

			glVertex3f(min_x, max_y, -0.1f); 
			glVertex3f(max_x, max_y, -0.1f); 

			glEnd( );
		}

	// draw circles
		glColor3f(1.0f, 1.0f, 1.0f);
		for (int i=0;i<number_of_circles;i++)
		{
			float cx = circles_center_x[i]; 
			float cy = circles_center_y[i]; 
			float r = circles_radius[i]; 

			glBegin(GL_LINES);
			for (int j = 0; j < 180; j++)
			{
			phi = j*pi/90;
			x = r * cos(phi);
			y = r * sin(phi);
			glVertex3f(cx + x, cy + y,0);
		    
			phi = (j+1)*pi/90;
			x = r * cos(phi);
			y = r * sin(phi);
			glVertex3f(cx + x, cy + y,0);
			}
			glEnd();
		}

	// draw rectangles
		glColor3f(1.0f, 1.0f, 1.0f);
		for (int i=0;i<number_of_rectangles;i++)
		{
			min_x = rectangles_min_x[i];
			min_y = rectangles_min_y[i];
			max_x = rectangles_max_x[i];
			max_y = rectangles_max_y[i];
			glBegin(GL_LINES);
			glVertex3f(min_x, min_y, -0.1f); 
			glVertex3f(min_x, max_y, -0.1f); 
			glVertex3f(min_x, min_y, -0.1f); 
			glVertex3f(max_x, min_y, -0.1f); 

			glVertex3f(max_x, min_y, -0.1f); 
			glVertex3f(max_x, max_y, -0.1f); 

			glVertex3f(min_x, max_y, -0.1f); 
			glVertex3f(max_x, max_y, -0.1f); 

			glEnd( );
		}

	// draw boundaries
			glColor3f(0.0f, 0.0f, 0.0f);
			min_x = domain_min_x;
			min_y = domain_min_y;
			max_x = domain_max_x;
			max_y = domain_max_y;
			glBegin(GL_LINES);
			glVertex3f(min_x, min_y, -0.1f); 
			glVertex3f(min_x, max_y, -0.1f); 
			glVertex3f(min_x, min_y, -0.1f); 
			glVertex3f(max_x, min_y, -0.1f); 

			glVertex3f(max_x, min_y, -0.1f); 
			glVertex3f(max_x, max_y, -0.1f); 

			glVertex3f(min_x, max_y, -0.1f); 
			glVertex3f(max_x, max_y, -0.1f); 

			glEnd( );

	// draw cpml boundaries
			glColor3f(1.0f, 1.0f, 0.0f);
			min_x = domain_min_x + n_cpml_xn*dx;
			min_y = domain_min_y + n_cpml_yn*dy;
			max_x = domain_max_x - n_cpml_xp*dx;
			max_y = domain_max_y - n_cpml_yp*dy;
			glBegin(GL_LINES);
			glVertex3f(min_x, min_y, -0.1f); 
			glVertex3f(min_x, max_y, -0.1f); 
			glVertex3f(min_x, min_y, -0.1f); 
			glVertex3f(max_x, min_y, -0.1f); 

			glVertex3f(max_x, min_y, -0.1f); 
			glVertex3f(max_x, max_y, -0.1f); 

			glVertex3f(min_x, max_y, -0.1f); 
			glVertex3f(max_x, max_y, -0.1f); 

			glEnd( );
		
	glEndList();
}

int getMaxDeviceId()
{
	int n_devices = 0;
	cudaGetDeviceCount(&n_devices);
	if (n_devices == 0) return false;
	int max_mem = 0;
	int max_flop = 0;
	int flop;
	int currentDevice = 0;
	int max_device = 0;
    cudaDeviceProp deviceProp;

	printf("<Available GPU devices on the system>\n", deviceProp.name);
	for (currentDevice=0;currentDevice<n_devices;currentDevice++)
	{
		cudaGetDeviceProperties(&deviceProp, currentDevice);
		flop = deviceProp.multiProcessorCount * deviceProp.clockRate;
		if (flop > max_flop)
		{
			max_device = currentDevice;
			max_flop = flop;
		} 
		printf("<%d> %s>\n", currentDevice, deviceProp.name);
		printf("  Processors: %d, Clock rate: %d, Global memory: %lu\n", 
			deviceProp.multiProcessorCount, deviceProp.clockRate,deviceProp.totalGlobalMem);
	}
	printf("<===================================>\n", deviceProp.name);
	return max_device;
}

bool runFdtdWithFieldDisplay(int argc, char** argv)
{
	// Initialize CUDA context
	cudaGLSetGLDevice (getMaxDeviceId());

	// Initialize GL context
	if( false == initializeGL(argc, argv)) return 0;

	// Initialize GL buffers
	initializeGLBuffers();

	// colormap used to map field intensity
	createColormapOnGpu();

	// Display list of objects in problem space
	createDisplayListForObjects();

	// copy data from CPU RAM to GPU global memory
	if (int ret = copyFdtdArraysToGpuMemory() !=0)
	{
		if (ret == 1) printf("Memory allocation error!\nProgram cannot run!\n");
		return 0;
	}

	time(&start);

	glutMainLoop();	// GLUT loop 

    Cleanup(0);
}

bool runFdtdWithoutFieldDisplayOnGpu()
{
	if (int ret = copyFdtdArraysToGpuMemory() !=0)
	{
		if (ret == 1) printf("Memory allocation error!\nProgram cannot run!\n");
		return 0;
	}

	time(&start);

	while (time_step<number_of_time_steps) fdtdIterationOnGpu();

	time(&stop);
	printf("Time elapsed in loop: %f sec\n ",(double )difftime(stop, start));


	fetchResultsFromGpuMemory();
	deallocateCudaArrays();
	deallocateArrays();
	saveSampledFieldsToFile();
}

void runFdtdWithoutFieldDisplayOnCpu()
{
	time(&start);

	while (time_step<number_of_time_steps) fdtdIterationOnCpu();

	time(&stop);
	printf("Time elapsed in loop: %f sec\n ",difftime(stop, start));

	deallocateArrays();
	saveSampledFieldsToFile();
}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
	if (!read_project_parameters_from_file(argc, argv)) return 0;
	if (!initialize_sources()) return 0;

	if (computation_platform == 2) // run on cpu
	{
		runFdtdWithoutFieldDisplayOnCpu();
	}
	if (computation_platform == 3) // run on gpu
	{
		if (show_Hz || show_Ez)
			runFdtdWithFieldDisplay(argc, argv);
		else
			runFdtdWithoutFieldDisplayOnGpu();
	}


}