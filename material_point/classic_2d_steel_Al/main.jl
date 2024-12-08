include("setup_works.jl")

using Printf
using LinearAlgebra: I, det
import PyPlot

include("material_point_type.jl")
include("mpm_grid.jl")
include("my_basis.jl")

#pyFig_RealTime = PyPlot.figure("MPM Disk impact", figsize=(16/2.54, 16/2.54), edgecolor="white", facecolor="white")

function mpmMain()
    fGravity = 0.0

    fOffset = 60.0/50/2.0
    thisMaterialDomain_01 = createMaterialDomain_Circle([30.0; 50.0], 9.6/2.0, fOffset)
    volume = fOffset*fOffset #3.14159*0.2*0.2/length(thisMaterialDomain_01)
    mass = 7850e-12*volume
    for thisMat in thisMaterialDomain_01
        thisMat.mass = mass
        thisMat.initial_volume = volume
        thisMat.volume = volume
        thisMat.elastic_modulus = 200.0e3
        thisMat.poisson_ratio = 0.3
        thisMat.yield_stress = 1.0e24
        thisMat.velocity = [0.0; -1160.0e3]
        thisMat.momentum = mass*thisMat.velocity
        thisMat.ext_force = [0.0; -fGravity*mass]
        thisMat.deform_grad[:,:] .= Matrix(I(2))
        thisMat.deform_grad_incr[:,:] .= Matrix(I(2))
    end

    thisMaterialDomain_02 = createMaterialDomain_Rectangle([30.0; 20.0], 60.0, 40.6, fOffset)
    volume = fOffset*fOffset #3.14159*0.2*0.2/length(thisMaterialDomain_02)
    mass = 2700e-12*volume
    for thisMat in thisMaterialDomain_02
        thisMat.mass = mass
        thisMat.initial_volume = volume
        thisMat.volume = volume
        thisMat.elastic_modulus = 78.2e3
        thisMat.poisson_ratio = 0.3
        thisMat.yield_stress = 300.0
        thisMat.velocity = [0.0; 0.0]
        thisMat.momentum = mass*thisMat.velocity
        thisMat.ext_force = [0.0; -fGravity*mass]
        thisMat.deform_grad[:,:] .= Matrix(I(2))
        thisMat.deform_grad_incr[:,:] .= Matrix(I(2))
    end

    Npoints1 = length(thisMaterialDomain_01)
    Npoints2 = length(thisMaterialDomain_02)
    NpointsAll = Npoints1 + Npoints2
    # array holding all material points (these are references to MaterialDomain_01 & 02)
    allMaterialPoint = Vector{mpmMaterialPoint_2D_Classic}(undef,NpointsAll)
    ip_all = 0
    for iIndex_MP in 1:Npoints1
        ip_all += 1
        allMaterialPoint[ip_all] = thisMaterialDomain_01[iIndex_MP]
    end
    for iIndex_MP in 1:Npoints2
        ip_all += 1
        allMaterialPoint[ip_all] = thisMaterialDomain_02[iIndex_MP]
    end

    # ---------------------------------------------------------------------------
    # information about the created domain
    mass = 0.0
    for iIndex_MP in 1:1:length(allMaterialPoint)
        mass += allMaterialPoint[iIndex_MP].mass
    end
    @printf("Initial configuration: \n")
    @printf("    Single particle Mass: %+.6e \n", allMaterialPoint[1].mass)
    @printf("    Total mass: %+.6e \n", mass)

    @printf("    Disk, number of material points: %d \n", length(thisMaterialDomain_01))
    @printf("    Target, number of material points: %d \n", length(thisMaterialDomain_02))
    @printf("    Total number of material points: %d \n", length(allMaterialPoint))


    # grid creation
    # nodes where fixation boundary conditions are hard coded in the following code!!!
    thisGrid = MPMGrid(60.0, 60.0, 51, 51)


    # ---------------------------------------------------------------------------
    # timers
    # ---------------------------------------------------------------------------
    # analysis timer
    fTimeIncrement = 1.0e-8
    fTimeEnd = 1.0e-4
    iTimeCycle = 0

    # realtime graphics timer
    fPlotTimeInterval = 100.0*fTimeIncrement
    fPlotTime = fPlotTimeInterval

    # final results plot timer
    fResultTimeInterval = fTimeIncrement#*fTimeEnd
    fResultTime = 0

    # console output timer
    fConsolTimeInterval = 0.1*fTimeEnd#1000.0*fTimeIncrement
    fConsolTime = fConsolTimeInterval

    # profiler timers
    fProfiler_MainLoop = 0.0
    fProfiler_Particle2Grid = 0.0
    fProfiler_Grid2Particle = 0.0
    # ---------------------------------------------------------------------------
    # plot arrays
    # ---------------------------------------------------------------------------
    # final results plot holder arrays
    plot_Time = Vector{Float64}(undef, 0)
    plot_Displacement = Array{Float64}(undef, 0)
    plot_KineticEnergy = Array{Float64}(undef, 0)
    plot_StrainEnergy = Array{Float64}(undef, 0)

    # main analysis loop
    #time_range = 0.0:fTimeIncrement:fTimeEnd
    time_range = 0.0:fTimeIncrement:2*fTimeIncrement
    @info "Length of time_range = $(length(time_range))"

    for fTime in time_range

        iTimeCycle += 1
        # ------------------------------------------------------------------------
        # realtime graphical plotting routines
        # @printf("Plotting...")
        # ------------------------------------------------------------------------
        fPlotTime += fTimeIncrement
        
        #=
        if(fPlotTime > fPlotTimeInterval)
            fPlotTime = 0.0

            iMaterialPoints = length(allMaterialPoint)
            array_x = [allMaterialPoint[i].centroid[1] for i in 1:iMaterialPoints]
            array_y = [allMaterialPoint[i].centroid[2] for i in 1:iMaterialPoints]
            array_color = Matrix{Float64}(undef, iMaterialPoints, 3)
            array_size = Matrix{Float64}(undef, iMaterialPoints, 1)
            
            # Set color ?
            for iIndex in 1:iMaterialPoints
                thisColor = allMaterialPoint[iIndex].equiv_plastic_strain
                thisColor /= (allMaterialPoint[iIndex].yield_stress/allMaterialPoint[iIndex].elastic_modulus)*500

                if(thisColor > 1.0)
                    # @printf("thiscolor %f\n", thisColor)
                    thisColor = 1.0
                end

                if(allMaterialPoint[iIndex].yield_stress > 1.0e3)
                    array_color[iIndex, :] = [1.0, 0.0, 0.0]
                else
                  array_color[iIndex, :] = [thisColor, 0.5*thisColor, 1.0-thisColor]#[thisGrid.points[iIndex].mass/iMaterialPoints, 0.0, 0.0]
                end
                array_size[iIndex, :] = [4.0]
           end

            pyPlot01 = PyPlot.gca()
            # pyPlot01 = PyPlot.subplot2grid((1,1), (0,0), colspan=1, rowspan=1, aspect="equal")
            PyPlot.scatter(array_x, array_y, c=array_color, lw=0, s=array_size)
            #
            pyPlot01[:spines]["top"][:set_color]("gray")
            pyPlot01[:spines]["right"][:set_color]("gray")
            pyPlot01[:spines]["bottom"][:set_color]("gray")
            pyPlot01[:spines]["left"][:set_color]("gray")
            # pyPlot01[:axhline](linewidth=4, color="g")
            # pyPlot01[:axvline](linewidth=4, color="g")
            pyPlot01[:set_xlim](0.0, 60.0)
            pyPlot01[:set_ylim](0.0, 60.0)
            # pyPlot01[:set_xlabel]("")
            # pyPlot01[:set_ylabel]("")
            pyPlot01[:grid](true, which="both", color="white", linestyle="-", linewidth=0.2)
            pyPlot01[:set_axisbelow](true)
            pyPlot01[:set_xticks]([])# empty to have no major ticks and grids
            pyPlot01[:set_xticks](collect(0.0:1.2:60.0),minor=true)
            pyPlot01[:set_yticks]([])# empty to have no major ticks and grids
            pyPlot01[:set_yticks](collect(0.0:1.2:60.0),minor=true)

            # PyPlot.show()
            # PyPlot.hold(true)

            strFileName = "TEMP_DiskImpact_$(iTimeCycle).png"
            PyPlot.savefig(strFileName, bbox_inches="tight")
            #PyPlot.hold(false)
        end
        =#

        #reset grid------------------------------------
        for iIndex in 1:thisGrid.NnodesTotal
            thisGrid.points[iIndex].mass = 0.0
            thisGrid.points[iIndex].velocity[:,:] .= 0.0
            thisGrid.points[iIndex].momentum[:,:] .= 0.0
            thisGrid.points[iIndex].force[:,:] .= 0.0
        end

        # material to grid -------------------------------------------------------
        for iIndex_MP in 1:NpointsAll
            thisMaterialPoint = allMaterialPoint[iIndex_MP]
            thisAdjacentGridPoints = getAdjacentGridPoints(allMaterialPoint[iIndex_MP], thisGrid)
            Nadjacent = length(thisAdjacentGridPoints)
            #
            for iIndex in 1:Nadjacent
                # sina, be careful here, this might not be by reference and might not be good for assignment
                thisGridPoint = thisGrid.points[thisAdjacentGridPoints[iIndex]]

                fShapeValue = getShapeValue_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

                v2ShapeGradient = getShapeGradient_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

                # mass
                thisGridPoint.mass += fShapeValue * thisMaterialPoint.mass
                # momentum
                thisGridPoint.momentum += fShapeValue * thisMaterialPoint.mass * thisMaterialPoint.velocity
                # internal forces
                volume = thisMaterialPoint.volume
                thisGridPoint.force[1] += -volume * (v2ShapeGradient[1]*thisMaterialPoint.stress[1] + v2ShapeGradient[2]*thisMaterialPoint.stress[3])
                thisGridPoint.force[2] += -volume * (v2ShapeGradient[2]*thisMaterialPoint.stress[2] + v2ShapeGradient[1]*thisMaterialPoint.stress[3])
                # external forces
                thisGridPoint.force += fShapeValue*thisMaterialPoint.ext_force
           end
        end
        #
        # update grid momentum and apply boundary conditions
        for iIndex_GP in 1:thisGrid.NnodesTotal
            thisGridPoint = thisGrid.points[iIndex_GP]
            thisGridPoint.momentum += thisGridPoint.force * fTimeIncrement
            if thisGridPoint.is_fixed[1]
                thisGridPoint.momentum[1] = 0.0
                thisGridPoint.force[1] = 0.0
            end
            if thisGridPoint.is_fixed[2]
                thisGridPoint.momentum[2] = 0.0
                thisGridPoint.force[2] = 0.0
            end
        end

        # ------------------------------------------------------------------------
        # grid to material pass 1-------------------------------------------------
        for iIndex_MP in 1:NpointsAll
            #
            thisMaterialPoint = allMaterialPoint[iIndex_MP]
            thisAdjacentGridPoints = getAdjacentGridPoints(thisMaterialPoint, thisGrid)
            #v2CentroidIncrement = zeros(2)
            for iIndex in 1:1:length(thisAdjacentGridPoints)
                thisGridPoint = thisGrid.points[thisAdjacentGridPoints[iIndex]]

                fShapeValue = getShapeValue_Classic(thisMaterialPoint, thisGridPoint, thisGrid)
                #v2ShapeGradient = getShapeGradient_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

                thisMaterialPoint.velocity += (fShapeValue * thisGridPoint.force / thisGridPoint.mass) * fTimeIncrement
           end
        end
        # reset grid momenta -----------------------------------------------------
        # for iIndex_GP in 1:1:thisGrid.NnodesTotal
        #     thisGrid.points[iIndex_GP].momentum = [0.0; 0.0]
        # end
        # map particle momenta back to grid------------------------------
        # mass in NOT mapped here
        for iIndex_MP in 1:1:length(allMaterialPoint)
            thisMaterialPoint = allMaterialPoint[iIndex_MP]
            thisAdjacentGridPoints = getAdjacentGridPoints(thisMaterialPoint, thisGrid)
            for iIndex in 1:1:length(thisAdjacentGridPoints)
                thisGridPoint = thisGrid.points[thisAdjacentGridPoints[iIndex]]

                fShapeValue = getShapeValue_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

                # thisGridPoint.momentum += fShapeValue * thisMaterialPoint.mass * thisMaterialPoint.velocity
                thisGridPoint.velocity += fShapeValue * thisMaterialPoint.mass * thisMaterialPoint.velocity / thisGridPoint.mass
            end
        end
        #apply boundary conditions velocity------------------------------------
        for iIndex_GP in 1:1:thisGrid.NnodesTotal
            thisGridPoint = thisGrid.points[iIndex_GP]

            if thisGridPoint.is_fixed[1]
                thisGridPoint.velocity[1] = 0.0
                thisGridPoint.momentum[1] = 0.0
                thisGridPoint.force[1] = 0.0
            end
            if thisGridPoint.is_fixed[2]
                thisGridPoint.velocity[2] = 0.0
                thisGridPoint.momentum[2] = 0.0
                thisGridPoint.force[2] = 0.0
            end
        end
        # ------------------------------------------------------------------------
        v2CentroidIncrement = zeros(Float64, 2)
        v3StrainIncrement = zeros(Float64, 3)

        # grid to material pass 2 ------------------------------------------------
        for iIndex_MP in 1:NpointsAll
            thisMaterialPoint = allMaterialPoint[iIndex_MP]
            thisAdjacentGridPoints = getAdjacentGridPoints(thisMaterialPoint, thisGrid)
            fill!(v2CentroidIncrement, 0.0)

            Nadjacent = length(thisAdjacentGridPoints)
            #println("Nadjacent = ", Nadjacent)
            for iIndex in 1:Nadjacent
                thisGridPoint = thisGrid.points[thisAdjacentGridPoints[iIndex]]

                fShapeValue = getShapeValue_Classic(thisMaterialPoint, thisGridPoint, thisGrid)
                v2ShapeGradient = getShapeGradient_Classic(thisMaterialPoint, thisGridPoint, thisGrid)

                v2GridPointVelocity = thisGridPoint.velocity#momentum / thisGridPoint.mass

                v2CentroidIncrement .+= (fShapeValue * thisGridPoint.momentum / thisGridPoint.mass) * fTimeIncrement

                # from (2011) A convected particle domain interpolation technique to extend ...
                thisMaterialPoint.deform_grad_incr += v2GridPointVelocity*transpose(v2ShapeGradient)*fTimeIncrement;
            end
            thisMaterialPoint.centroid += v2CentroidIncrement
            thisMaterialPoint.deform_grad = thisMaterialPoint.deform_grad_incr * thisMaterialPoint.deform_grad
            v3StrainIncrement[1] = thisMaterialPoint.deform_grad_incr[1,1] - 1.0
            v3StrainIncrement[2] = thisMaterialPoint.deform_grad_incr[2,2] - 1.0
            v3StrainIncrement[3] = thisMaterialPoint.deform_grad_incr[1,2] + thisMaterialPoint.deform_grad_incr[2,1]
            thisMaterialPoint.deform_grad_incr[:,:] .= Matrix(I(2))
            thisMaterialPoint.strain[1] += v3StrainIncrement[1]
            thisMaterialPoint.strain[2] += v3StrainIncrement[2]
            thisMaterialPoint.strain[3] += v3StrainIncrement[3]

            fE = thisMaterialPoint.elastic_modulus;
            fNu = thisMaterialPoint.poisson_ratio
            fYield = thisMaterialPoint.yield_stress

            v3StrainCurrent = thisMaterialPoint.strain
            v3StressCurrent = thisMaterialPoint.stress
            v3PlasticStrainCurrent = thisMaterialPoint.plastic_strain
            fAlphaCurrent = thisMaterialPoint.equiv_plastic_strain

            v32Result = getIncrement_Plastic(fE, fNu, fYield, fAlphaCurrent, v3StressCurrent, v3StrainCurrent, v3PlasticStrainCurrent, v3StrainIncrement)

            v3StressIncrement = v32Result[:, 1]
            v3PlasticStrainIncrement = v32Result[:, 2]
            fAlphaIncrement = v32Result[1, 3]

            thisMaterialPoint.stress += v3StressIncrement
            thisMaterialPoint.plastic_strain += v3PlasticStrainIncrement
            thisMaterialPoint.equiv_plastic_strain += fAlphaIncrement

            thisMaterialPoint.volume = det(thisMaterialPoint.deform_grad) * thisMaterialPoint.initial_volume

            thisMaterialPoint.momentum = thisMaterialPoint.velocity * thisMaterialPoint.mass
        end

        # ------------------------------------------------------------------------
        # calculating strain and kinetic energy for final results plot
        # ------------------------------------------------------------------------
        fResultTime += fTimeIncrement
        if fResultTime > fResultTimeInterval
            fResultTime = 0.0
            fStrainEnergy = 0.0
            fKineticEnergy = 0.0
            for iIndex_MP in 1:NpointsAll
                thisMaterialPoint = allMaterialPoint[iIndex_MP]
                fStrainEnergy += 0.5*thisMaterialPoint.strain[1] * thisMaterialPoint.stress[1] * thisMaterialPoint.volume
                fStrainEnergy += 0.5*thisMaterialPoint.strain[2] * thisMaterialPoint.stress[2] * thisMaterialPoint.volume
                fStrainEnergy += 0.5*thisMaterialPoint.strain[3] * thisMaterialPoint.stress[3] * thisMaterialPoint.volume

                fVelocity = thisMaterialPoint.velocity[1]^2 + thisMaterialPoint.velocity[2]^2
                fKineticEnergy += 0.5*fVelocity * thisMaterialPoint.mass
            end
            @info "$(fTime) $(fKineticEnergy) $(fStrainEnergy)"
            #save to plot arrays
            #push!(plot_Time, fTime)
            #push!(plot_KineticEnergy, fKineticEnergy)
            #push!(plot_StrainEnergy, fStrainEnergy)
        end
        # ------------------------------------------------------------------------
        # consol output
        # ------------------------------------------------------------------------
        fConsolTime += fTimeIncrement
        if(fConsolTime > fConsolTimeInterval)
            fConsolTime = 0.0

            mass = 0.0
            fMomentum_x = 0.0
            for iIndex_MP in 1:Npoints1
                mass += thisMaterialDomain_01[iIndex_MP].mass
                fMomentum_x += thisMaterialDomain_01[iIndex_MP].momentum[1]
            end
            
            #fProfiler_Total = fProfiler_Particle2Grid + fProfiler_Grid2Particle
            #@printf("fTime: %+.3e |", fTime)
            #@printf("M_x: %+.3e |", fMomentum_x)
            #@printf("(Profiler) Total: %+.3e ", fProfiler_Total)
            #@printf("P2G: %+.3e (%+.2f) ", fProfiler_Particle2Grid, fProfiler_Particle2Grid/fProfiler_Total)
            #@printf("G2P: %+.3e (%+.2f) \n", fProfiler_Grid2Particle, fProfiler_Grid2Particle/fProfiler_Total)
        end
    end

    # final plots
    # pyFig_RealTime = PyPlot.figure("2Disk Impact FinalPlot", figsize=(8/2.54, 4/2.54))
    # PyPlot.clf()
    # pyPlot01 = PyPlot.gca()
    # PyPlot.subplots_adjust(left=0.15, bottom=0.25, right=0.65)
    # pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
    # pyPlot01[:set_axisbelow](true)
    # pyPlot01[:set_xlim](0.0, 4.0)
    # pyPlot01[:set_ylim](0.0, 3.0)
    # pyPlot01[:set_xlabel]("time (s)", fontsize=8)
    # pyPlot01[:set_ylabel]("energy (\$\\times 10^{-3}\$ Nm)", fontsize=8)
    # pyPlot01[:set_xticks](collect(0.0:1.0:4.0))
    # pyPlot01[:tick_params](axis="both", which="major", labelsize=8)
    # pyPlot01[:set_yticks](collect(0.0:1.0:3.0))
    # PyPlot.plot(plot_Time, c="blue", plot_KineticEnergy, "-", label="\$ K \$", linewidth=1.0)
    # PyPlot.hold(true)
    # PyPlot.plot(plot_Time, c="red", plot_StrainEnergy, "-", label="\$ U \$", linewidth=1.0)
    # PyPlot.plot(plot_Time, c="green", plot_KineticEnergy + plot_StrainEnergy, "-", label="\$ K+U \$", linewidth=1.0)
    # PyPlot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=8)
    # PyPlot.savefig("..\\..\\Figs\\plot_2Disk_Julia.pdf")

#=
=#


end

@time mpmMain()
