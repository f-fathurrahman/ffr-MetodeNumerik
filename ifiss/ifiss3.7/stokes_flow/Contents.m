% STOKES_FLOW
%
% Files
%   altinfsup               - inf-sup eigenvalue distribution of potential flow
%   box_stokes              - set up flow problem in box domain
%   channel_stokes          - set up inflow/outflow problem in channel domain
%   flowbc                  - imposes inflow boundary condition
%   forwardstep_stokes      - set up flow problem in forward step domain
%   helpme_stokes           - Stokes flow problem interactive help
%   helpme_twofieldpressure - Stokes flow problem interactive help
%   infsup                  - computes inf-sup eigenvalue distribution
%   localbc_xy              - imposes vector BC for Poisson error estimator
%   longstep_stokes         - set up flow problem in extended step domain
%   obstacle_stokes         - set up flow problem in domain with square obstacle
%   out_Neumann_bound       - locates elements on the step outflow boundary
%   plate_stokes            - set up flow problem in slit domain
%   q1div                   - computes norm of divergence of Q1 flow solution
%   q1divz                  - legacy code | replaced by corrected q1div.m
%   q2div                   - computes norm of divergence of Q2 flow solution
%   q2divz                  - legacy code | replaced by corrected q2div.m
%   quad_stokes             - set up Stokes problem in quadrilateral domain
%   solve_stokes            - solve Stokes problem in square domain
%   solve_step_stokes       - solve Stokes problem in step domain
%   solveX_stokes           - solve two-field Stokes problem in square domain
%   solveX_step_stokes      - solve two-field Stokes problem in step domain
%   specific_flow           - (current) problem flow boundary condition
%   square_stokes           - set up flow problem in unit square domain
%   squareX_stokes          - set up two-field approximation in unit square domain
%   step_stokes             - set up flow problem in standard step domain
%   stepX_stokes            - set up two-field approximation in step domain
%   stokes_q1p0             - Q1-P0 matrix generator
%   stokes_q1q1             - Q1-Q1 matrix generator
%   stokes_q2p1             - Q2-P1 matrix generator
%   stokes_q2q1             - Q2-Q1 matrix generator
%   stokes_q2q1q0           - Q2-Q1* matrix generator
%   stokespost              - estimates Stokes error distribution
%   stokespost_q1p0_bc      - postprocesses Poisson error estimator
%   stokespost_q1p0_p       - computes Poisson error estimator for Q1-P0
%   stokespost_q1q1_bc      - postprocesses Poisson error estimator
%   stokespost_q1q1_p       - computes Poisson error estimator for Q1-Q1
%   stokespost_q2p1         - vectorized Q2-P1 local Poisson error estimator
%   stokespost_q2q1         - vectorized Q2-Q1 local Poisson error estimator
%   stream_bc               - (current) streamfunction boundary condition
%   streambc                - imposes Dirichlet BC on the streamfunction
%   stressjmps_q1p0         - stress jumps for rectangular Q1-P0 grid
%   stressjmps_q1q1         - stress jumps for rectangular Q1-Q1 grid
%   stressjmps_q2p1         - vectorised stress jumps for Q2-P1 grid
%   stressjmps_q2q1         - vectorised stress jumps for Q2-Q1 grid
%   symstep_stokes          - set up flow problem in symmetric step domain
%   vorticity_q1            - Q1 total vorticity
%   vorticity_q2            - Q2 total vorticity
%   xsolve_stokes           - solve Stokes problem in channel domain
