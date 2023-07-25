!#define _TIMER_
!#define _HPC_FUGAKU_
!#define _READ_FROM_BUCKET_
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose the problem dimension
! 2D > _2D_ TODO
! 3D > _3D_
#define _3D_
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose the flow configuration
! Homogeneous triperiodic flow  > _FLOW_TRIPERIODIC_
! Plane Couette flow            > _FLOW_COUETTE_
! Plane Poiseuille channel flow > _FLOW_CHANNEL_
! Fixed flow                    > _FLOW_FIXED_
! Flow with Inlet-Outlet        > _FLOW_INOUT_
#define _FLOW_TRIPERIODIC_

#ifdef _FLOW_TRIPERIODIC_
! For _FLOW_TRIPERIODIC_ choose the driving force
! Taylor-Green vortex flow        > _FLOW_TRIPERIODIC_TGV_
! Arnold-Beltrami-Childress flow  > _FLOW_TRIPERIODIC_ABC_
! Kolmogorov parallel flow        > _FLOW_TRIPERIODIC_KOL_
! Homogeneous isotropic forcing   > _FLOW_TRIPERIODIC_ISO_
! Linear forcing in phys. space   > _FLOW_TRIPERIODIC_LIN_
#define _FLOW_TRIPERIODIC_ISO_
#ifdef _FLOW_TRIPERIODIC_KOL_
! Kolmogorov-like parallel flow with mean velocity refVel in (choose one):
! x direction > _FLOW_TRIPERIODIC_KOL_MEANU_
! y direction > _FLOW_TRIPERIODIC_KOL_MEANV_
! z direction > _FLOW_TRIPERIODIC_KOL_MEANW_
!#define _FLOW_TRIPERIODIC_KOL_MEANV_
#endif
! write mean fluid velocity  > _WRITE_MEAN_VELOCITY_
!#define _WRITE_MEAN_VELOCITY_
! remove mean fluid velocity > _REMOVE_MEAN_VELOCITY_
#define _REMOVE_MEAN_VELOCITY_
#endif

#ifdef _FLOW_CHANNEL_
! For _FLOW_CHANNEL_ choose the driving force
! Constant flow rate         > _FLOW_CHANNEL_CFR_
! Constant pressure gradient > _FLOW_CHANNEL_CPG_
#define _FLOW_CHANNEL_CFR_
#ifdef _FLOW_CHANNEL_CPG_
! For _FLOW_CHANNEL_CPG_ choose the shape of the driving force
! Constant pressure gradient  > _FLOW_CHANNEL_CPG_STD_
! Pulsating pressure gradient > _FLOW_CHANNEL_CPG_PLS_
#define _FLOW_CHANNEL_CPG_PLS_
#endif
! For _FLOW_CHANNEL_ if simulating only half of channel
! Half channel > _FLOW_CHANNEL_HALF_
! Full channel > _FLOW_CHANNEL_FULL_
#define _FLOW_CHANNEL_FULL_
! For _FLOW_CHANNEL_ activate array of probes
! Flow probes > _FLOW_CHANNEL_PROB_
! No probes   > _FLOW_CHANNEL_NOPR_
!#define _FLOW_CHANNEL_NOPR_
#endif

#if defined(_FLOW_TRIPERIODIC_) || defined(_FLOW_CHANNEL_)
! For _FLOW_TRIPERIODIC_ and _FLOW_CHANNEL_ choose the initial condtion
! At rest   > _FLOW_INIT_REST_
! Laminar   > _FLOW_INIT_LAMINAR_
! Turbulent > _FLOW_INIT_TURBULENT_
! Wave (only when vof and interface) > _POTENTIAL_WAVE_
#define _FLOW_INIT_LAMINAR_
#endif

#ifdef _FLOW_INOUT_
! For _FLOW_INOUT_ choose the flow configuration
! Channel flow                           > _FLOW_INOUT_CHANNEL_
! Unbounded open flow                    > _FLOW_INOUT_UNBOUND_
! Unbounded open flow (periodic along Z) > _FLOW_INOUT_UNBOUND_PER_
! Jet                                    > _FLOW_INOUT_JETPUFF_
#define _FLOW_INOUT_JETPUFF_
! define position of probes
! probes only at centerline              > _FLOW_INOUT_PRB_CEN_
! probes also around centerline          > _FLOW_INOUT_PRB_LAT_
#define _FLOW_INOUT_PRB_CEN_
#ifdef _FLOW_INOUT_JETPUFF_
! For _FLOW_INOUT_JETPUFF_ choose planar or round jet
! planar jet          > _FLOW_INOUT_JETPUFF_PLN_
! round jet           > _FLOW_INOUT_JETPUFF_RND_
#define _FLOW_INOUT_JETPUFF_PLN_
! For _FLOW_INOUT_ choose whether inlet condition is constant over time or uses time-dependent inlet condition
! constant inlet        > (default)
! time-dependent inlet  > _FLOW_INOUT_PUFF_TIME_
!#define _FLOW_INOUT_PUFF_TIME_
! For _FLOW_INOUT_JETPUFF_ choose initial condition
! zero-flow             > _FLOW_INOUT_JETPUFF_INIZER_
! plug flow (cylinder)  > _FLOW_INOUT_JETPUFF_INICYL_
! sech^2 flow           > _FLOW_INOUT_JETPUFF_INISEC_
#define _FLOW_INOUT_JETPUFF_INISEC_
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose the pressure solver
! Fourier in X and Y           > _SOLVER_FFT2D_
! Fourier in X, Y and Z        > _SOLVER_FFT3D_
! Fourier in X and Cosine in Y > _SOLVER_INOUT_
#define _SOLVER_FFT3D_
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _SOLVER_FFT2D_
! Add stretching in Z
! Stretched grid implemented only for _SINGLEPHASE_, _MULTIPHASE_VOF_ (can't use with _MULTIPHASE_SURFACTANT_), _MULTIPHASE_IBML_
! Z uniform grid   > _SOLVER_FFT2D_UNIZ_
! Z stretched grid > _SOLVER_FFT2D_STRZ_
#define _SOLVER_FFT2D_UNIZ_
! Be careful: missing correct bulkVel for _FLOW_TRIPERIODIC_KOL_ if you want to use it with FFT2D and Z stretched
#ifdef _SOLVER_FFT2D_STRZ_
! Define the stretched Z grid
! Read from file > _SOLVER_FFT2D_STRZ_READ_
! Compute (tanh) > _SOLVER_FFT2D_STRZ_COMP_
#define _SOLVER_FFT2D_STRZ_COMP_
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose the multiphase system
! Singlephase (uniform viscosity and density) > _SINGLEPHASE_
! Multiphase                                  > _MULTIPHASE_
#define _MULTIPHASE_

#ifdef _MULTIPHASE_
! For _MULTIPHASE_ choose the kind of flow
! Non-Newtonian  > _MULTIPHASE_NNF_
! Particles      > _MULTIPHASE_IBM_ !XXX change name
! Lagrangian IBM > _MULTIPHASE_IBML_
! Droplets       > _MULTIPHASE_VOF_
! Surfactant     > _MULTIPHASE_SURFACTANT_ (flag inside ifdef(VoF), no VoF, no surfactant)
#define _MULTIPHASE_VOF_
! For _MULTIPHASE_ choose the viscosity treatment !IC only for _MULTIPHASE_VOF_ right?
! Uniform viscosity  > _MULTIPHASE_VISCOSITY_UNIFORM_
! Variable viscosity > _MULTIPHASE_VISCOSITY_VARIABLE_
#define _MULTIPHASE_VISCOSITY_UNIFORM_
! For _MULTIPHASE_ choose the density treatment
! Uniform density  > _MULTIPHASE_DENSITY_UNIFORM_
! Variable density > _MULTIPHASE_DENSITY_VARIABLE_
#define _MULTIPHASE_DENSITY_UNIFORM_

#ifdef _MULTIPHASE_NNF_
! For _MULTIPHASE_NNF_ choose the fluid type
! Power-law fluid (Carreau) > _MULTIPHASE_NNF_POW_
! Oldroyd-B model           > _MULTIPHASE_NNF_OLB_
! PTT model                 > _MULTIPHASE_NNF_PTT_
! Giesekus model            > _MULTIPHASE_NNF_GIE_
! FENE-P model              > _MULTIPHASE_NNF_FEN_
! EVP Saramito              > _MULTIPHASE_NNF_EVP_
#define _MULTIPHASE_NNF_OLB_

#ifdef _MULTIPHASE_NNF_POW_
! For _MULTIPHASE_NNF_POW_ choose between thinning and thickening behaviour
! Shear-thinning   > _MULTIPHASE_NNF_POW_N_
! Shear-thickening > _MULTIPHASE_NNF_POW_C_
#define _MULTIPHASE_NNF_POW_C_
#endif

#ifdef _MULTIPHASE_NNF_EVP_
! EVP Saramito                   > _MULTIPHASE_NNF_EVP_SRM_
! EVP Saramito/PTT               > _MULTIPHASE_NNF_EVP_SRM_PTT_
! EVP Saramito/Herschel-Buckley  > _MULTIPHASE_NNF_EVP_SRM_HB_
#define _MULTIPHASE_NNF_EVP_SRM_HB_
#endif

#ifndef _MULTIPHASE_NNF_POW_
! For _MULTIPHASE_NNF_XXX_ choose between the stress or conformation tensor formulation
! Stress tensor       > _MULTIPHASE_NNF_XXX_TAU_
! Conformation tensor > _MULTIPHASE_NNF_XXX_CON_
! LogConformation     > _MULTIPHASE_NNF_XXX_LOG_
#define _MULTIPHASE_NNF_XXX_LOG_
#endif
#ifdef _MULTIPHASE_VOF_
! if VOF is active, decide which phase is non-Newtonian
! phi=0 (carrier) non-Newtonian   > _MULTIPHASE_VOF_NN_CAR_
! phi=1 (drop/jet) non-Newtonian  > _MULTIPHASE_VOF_NN_DRP_
#define _MULTIPHASE_VOF_NN_CAR_
#endif
#endif

#ifdef _MULTIPHASE_VOF_
! For _MULTIPHASE_VOF_ choose the problem type
! droplets             > _MULTIPHASE_VOF_DROPLETS_
! interface            > _MULTIPHASE_VOF_INTERFACE_
! Zalesak slotted disk > _MULTIPHASE_VOF_ZALESAK_
#define _MULTIPHASE_VOF_DROPLETS_
#ifdef _MULTIPHASE_VOF_DROPLETS_
! For _MULTIPHASE_VOF_DROPLETS_ write stats_droplet_XXX.out, which is only valid for single drop
! single drop          > _MULTIPHASE_VOF_DROPLETS_SINGLE_
#define _MULTIPHASE_VOF_DROPLETS_SINGLE_
#endif
#ifdef _MULTIPHASE_VOF_INTERFACE_
! For interface that changes position in the x, y direction > _INTERFACE_XY_DEPENDENT_
#define _INTERFACE_XY_DEPENDENT_
#endif
! For _MULTIPHASE_VOF_ choose the problem dimension
! 2D > _MULTIPHASE_VOF_2D_
! 3D > _MULTIPHASE_VOF_3D_
! VOF+NNT works only with VOF3D
#ifndef _MULTIPHASE_NNF_
!#define _MULTIPHASE_VOF_2D_
#endif
#define _MULTIPHASE_VOF_3D_
! For _MULTIPHASE_VOF_ choose the order of the interface
! linear    > _MULTIPHASE_VOF_LINEAR_
! quadratic > _MULTIPHASE_VOF_QUADRATIC_
#define _MULTIPHASE_VOF_QUADRATIC_

!#define _MULTIPHASE_SURFACTANT_
#ifdef _MULTIPHASE_SURFACTANT_
! type of surfactant: ... to be done
! insoluble              > _MULTIPHASE_SURF_INS_
! soluble in phase 0     > _MULTIPHASE_SURF_S0_
! soluble in phase 1     > _MULTIPHASE_SURF_S1_
! soluble in both phases > _MULTIPHASE_SURF_S01_
#define _MULTIPHASE_SURF_INS_
! surfactant initial condition for interfacial surfactant concentration
! linear z-direction surfactant concentration > _MULTIPHASE_SURFACTANT_INIT_LINEAR_
! diffusion test, interfacial                 > _MULTIPHASE_SURFACTANT_INIT_DIFFUSION_
! constant concentration                      > _MULTIPHASE_SURFACTANT_INIT_CONSTANT_
#define _MULTIPHASE_SURFACTANT_INIT_LINEAR_
#endif

!#define _MULTIPHASE_SURFACTANT_PF_
#ifdef _MULTIPHASE_SURFACTANT_PF_
! choose initial condition, default is the uniform surfactant concentration
!#define _MULTIPHASE_SURFACTANT_PF_INI_MXINT_
#define _MULTIPHASE_SURFACTANT_PF_INI_EQ_
!#define _MULTIPHASE_SURFACTANT_PF_INI_DIFFTEST_
#endif

#endif

#ifdef _MULTIPHASE_IBM_
! For using the binning in the collosion model > _COLLISION_BIN_
#define _COLLISION_BIN_
! For paralelising the particle time advancement > _PARALLEL_PARTICLES_
#define _PARALLEL_PARTICLES_
#endif


#ifdef _MULTIPHASE_IBM_
! For _MULTIPHASE_IBM_ choose the problem dimension
! 2D > _MULTIPHASE_IBM_2D_
! 3D > _MULTIPHASE_IBM_3D_
#define _MULTIPHASE_IBM_3D_
!For static objects > _MULTIPHASE_IBM_STATIC_
!#define _MULTIPHASE_IBM_STATIC_
#ifdef _MULTIPHASE_IBM_STATIC_
! For static particles (spheres or cylinders) > _MULTIPHASE_IBM_STATIC_PARTICLES_
! For static user-defined static shapes (initialized in init.f90) > _MULTIPHASE_IBM_STATIC_UD_SHAPE_
#define _MULTIPHASE_IBM_STATIC_UD_SHAPE_
#endif
! Integrate particle collisions with Crank-Nicholson > _MULTIPHASE_IBM_COL_CRANK_NIC_
#define _MULTIPHASE_IBM_COL_CRANK_NIC_
#endif

#ifdef _MULTIPHASE_IBML_
! For _MULTIPHASE_IBML_ choose the method
! RKPM      > _MULTIPHASE_IBML_RKPM_
! Goldstein > _MULTIPHASE_IBML_GLDS_
#define _MULTIPHASE_IBML_GLDS_
#ifdef _MULTIPHASE_IBML_RKPM_
! For _MULTIPHASE_IBML_RKPM_ choose the method
! Look-up table > _MULTIPHASE_IBML_RKPM_LUT_
! Standard      > _MULTIPHASE_IBML_RKPM_STD_
#define _MULTIPHASE_IBML_RKPM_STD_
#endif
#ifdef _MULTIPHASE_IBML_GLDS_
! 2D (Goldstein only!) > _MULTIPHASE_IBML_GLDS_2D_
! 3D (Goldstein only!) > _MULTIPHASE_IBML_GLDS_3D_
#define _MULTIPHASE_IBML_GLDS_3D_
#endif
! For _MULTIPHASE_IBML_ choose how to obtain the initial objects positions
! read from file                            > _MULTIPHASE_IBML_XXXX_READ_
! compute it somehow (not implemented yet!) > _MULTIPHASE_IBML_XXXX_COMP_
#define _MULTIPHASE_IBML_XXXX_COMP_
! For _MULTIPHASE_IBML_ choose to simulate a filament or any other non-moving object
! filament           > _MULTIPHASE_IBML_XXXX_FIL_
! rigid flap in xz   > _MULTIPHASE_IBML_XXXX_FLAPXZ_
! other fixed object > _MULTIPHASE_IBML_XXXX_OFO_
#define _MULTIPHASE_IBML_XXXX_FIL_
#ifdef _MULTIPHASE_IBML_XXXX_FLAPXZ_
! For _MULTIPHASE_IBML_XXXX_FLAPXZ_ choose to apply an active moment or not
! active  > _MULTIPHASE_IBML_XXXX_FLAPXZ_ACT_
! passive > _MULTIPHASE_IBML_XXXX_FLAPXZ_PAS_
#define _MULTIPHASE_IBML_XXXX_FLAPXZ_ACT_
#endif
#ifdef _MULTIPHASE_IBML_XXXX_FIL_
#ifdef _MULTIPHASE_IBML_RKPM_
! For _MULTIPHASE_IBML_XXXX_FIL_ choose the problem dimension
! 2D (RKPM only!) > _MULTIPHASE_IBML_RKPM_FIL_2D_
! 3D (RKPM only!) > _MULTIPHASE_IBML_RKPM_FIL_3D_
#define _MULTIPHASE_IBML_RKPM_FIL_3D_
#endif
! For _MULTIPHASE_IBML_XXXX_FIL_ choose the filaments boundary conditions
! free    > _MULTIPHASE_IBML_XXXX_FIL_FREE_
! hinged  > _MULTIPHASE_IBML_XXXX_FIL_HING_
! clamped > _MULTIPHASE_IBML_XXXX_FIL_CLAM_
#define _MULTIPHASE_IBML_XXXX_FIL_FREE_
#ifndef _MULTIPHASE_IBML_XXXX_READ_
! For _MULTIPHASE_IBML_XXXX_FIL_ choose how to obtain the initial filament positions
! fixed arrangement > _MULTIPHASE_IBML_XXXX_FIL_FXAR_
! generate random   > _MULTIPHASE_IBML_XXXX_FIL_RAND_
#define _MULTIPHASE_IBML_XXXX_FIL_RAND_
#endif
! For _MULTIPHASE_IBML_XXXX_FIL_ choose the filament-filament interaction model
! minimal collision model > _MULTIPHASE_IBML_XXXX_FIL_COLY_
! no collision model      > _MULTIPHASE_IBML_XXXX_FIL_COLN_
#define _MULTIPHASE_IBML_XXXX_FIL_COLY_
#ifdef _MULTIPHASE_IBML_XXXX_FIL_COLY_
! For _MULTIPHASE_IBML_XXXX_FIL_COLY_ choose about counting the collisions
! no counter                                     > _MULTIPHASE_IBML_XXXX_FIL_COLY_NOCOUNT_
! activate counter for fiber-to-fiber collisions > _MULTIPHASE_IBML_XXXX_FIL_COLY_COUNT_
#define _MULTIPHASE_IBML_XXXX_FIL_COLY_COUNT_
#endif
#endif
#endif
#if defined(_MULTIPHASE_IBML_XXXX_OFO_) && defined(_MULTIPHASE_IBML_GLDS_3D_)
! For _MULTIPHASE_IBML_GLDS_OFO_ choose dimension of Lagrangian elements
! 1-D (segments)  > _MULTIPHASE_IBML_GLDS_OFO_LIN_
! 2-D (triangles) > _MULTIPHASE_IBML_GLDS_OFO_SUR_
#define _MULTIPHASE_IBML_GLDS_OFO_LIN_
! Save 3D data files containing the forces due to the multiphase parts > _VERBOSE_
!#define _VERBOSE_
#endif
#endif

! For runtime Lagrangian tracers and FTLE postprocessing
! Yes Lagrangian postprocessing > _POSTP_CFD2LCSY_
! No Lagrangian postprocessing > _POSTP_CFD2LCSN_
!#define _POSTP_CFD2LCSY_
#define _POSTP_CFD2LCSN_

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
