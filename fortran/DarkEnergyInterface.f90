    module DarkEnergyInterface
    use precision
    use interpolation
    use classes
    implicit none

    private

    type, extends(TCambComponent) :: TDarkEnergyModel
        logical :: is_cosmological_constant = .false.
        integer :: num_perturb_equations = 0
        real(dl) :: wwmg = -1
        real(dl) :: wawmg = 0
        
    contains
    procedure :: Init
    procedure :: BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy !Get density perturbation and heat flux for sources
    procedure :: diff_rhopi_Add_Term
    procedure :: PerturbationInitial
    procedure :: PerturbationEvolve
    procedure :: PrintFeedback
    ! do not have to implement w_de or grho_de if BackgroundDensityAndPressure is inherited directly
    procedure :: w_de
    procedure :: grho_de
    procedure :: Effective_w_wa !Used as approximate values for non-linear corrections
    end type TDarkEnergyModel

    type, extends(TDarkEnergyModel) :: TDarkEnergyEqnOfState
        !Type supporting w, wa or general w(z) table
        real(dl) :: w_lam = -1_dl !p/rho for the dark energy (an effective value, used e.g. for halofit)
        real(dl) :: wa = 0._dl !may not be used, just for compatibility with e.g. halofit
        real(dl) :: cs2_lam = 1_dl !rest-frame sound speed, though may not be used
        logical :: use_tabulated_w = .false.  !Use interpolated table; note this is quite slow.
        logical :: no_perturbations = .false. !Don't change this, no perturbations is unphysical
        !Interpolations if use_tabulated_w=.true.
        Type(TCubicSpline) :: equation_of_state, logdensity
    contains
    procedure :: ReadParams => TDarkEnergyEqnOfState_ReadParams
    procedure :: Init => TDarkEnergyEqnOfState_Init
    procedure :: SetwTable => TDarkEnergyEqnOfState_SetwTable
    procedure :: PrintFeedback => TDarkEnergyEqnOfState_PrintFeedback
    end type TDarkEnergyEqnOfState

    public TDarkEnergyModel, TDarkEnergyEqnOfState
    contains

    function w_de(this, a,M,N,omega_de,omega_node,wM,cri,E2)
    class(TDarkEnergyModel) :: this
    real(dl) :: w_de, al
    real(dl), intent(IN) :: a,M,N,omega_de,omega_node,wM,cri,E2

!     w_de = omega_node*(1._dl+wM)*M/(omega_de*sqrt(omega_node*E2+N)) - 1._dl


    end function w_de  ! equation of state of the PPF DE

    function grho_de(this, a)  !relative density (8 pi G a^4 rho_de /grhov)
    class(TDarkEnergyModel) :: this
    real(dl) :: grho_de, al, fint
    real(dl), intent(IN) :: a
    
    end function grho_de

    subroutine PrintFeedback(this, FeedbackLevel)
    class(TDarkEnergyModel) :: this
    integer, intent(in) :: FeedbackLevel

    end subroutine PrintFeedback


    subroutine Init(this, State)
    use classes
    class(TDarkEnergyModel), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    end subroutine Init

    subroutine BackgroundDensityAndPressure(this, grhov, a, grhov_t, w,M,N,grho_node_t,gpre_node_t,cri)
    !Get grhov_t = 8*pi*rho_de*a**2 and equation of state at scale factor a
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(in) :: grhov,a,  M,N,grho_node_t,cri,gpre_node_t
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w
    real(dl) :: omega_node,omega_de,A1,grhoN,grhoM2,grhoA1,wM,E2


!!!!WMG!!!!

    grhoM2 = M**2._dl*cri
    grhoN = N*cri
      
    grhov_t= (grhov - 2._dl*M*sqrt(cri)*sqrt(cri-grhov+grhoN) ) *a**2._dl + 2._dl*M*sqrt(cri)* sqrt(grho_node_t+(grhoN)*a**4._dl)

!!!!WMG!!!!    

!     omega_node= grho_node_t/cri/a**4
!     omega_de = grhov_t/cri/a**2
! 
!     wM = gpre_node_t/grho_node_t
!     
!     
!     E2= (grhov_t/a**2+ grho_node_t/a**4)/cri
!     w = this%w_de(a,M,N,omega_de,omega_node,wM,cri,E2)
!     this%wwmg=this%w_de(1._dl,M,N,omega_de,omega_node,wM,cri,1._dl)
!     this%wawmg=-M
    end subroutine BackgroundDensityAndPressure

    subroutine Effective_w_wa(this, w, wa)
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(out) :: w, wa

    w = this%wwmg
    wa = this%wawmg

    end subroutine Effective_w_wa


    subroutine PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe=0
    dgqe=0

    end subroutine PerturbedStressEnergy


    function diff_rhopi_Add_Term(this, dgrhoe, dgqe,grho, gpres, w, grhok, adotoa, &
        Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    class(TDarkEnergyModel), intent(in) :: this
    real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, grhok, w, adotoa, &
        k, grhov_t, z, k2, yprime(:), y(:), Kf1
    integer, intent(in) :: w_ix
    real(dl) :: ppiedot

    ! Ensure, that the result is set, when the function is not implemented by
    ! subclasses
    ppiedot = 0._dl

    end function diff_rhopi_Add_Term

    subroutine PerturbationEvolve(this, ayprime, w, w_ix, a, adotoa, k, z, y)
    class(TDarkEnergyModel), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a,adotoa, k, z, y(:), w
    integer, intent(in) :: w_ix
    end subroutine PerturbationEvolve

    subroutine PerturbationInitial(this, y, a, tau, k)
    class(TDarkEnergyModel), intent(in) :: this
    real(dl), intent(out) :: y(:)
    real(dl), intent(in) :: a, tau, k
    !Get intinitial values for perturbations at a (or tau)
    !For standard adiabatic perturbations can usually just set to zero to good accuracy

    y = 0

    end subroutine PerturbationInitial


    subroutine TDarkEnergyEqnOfState_SetwTable(this, a, w, n)
    class(TDarkEnergyEqnOfState) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n), w(n)
    real(dl), allocatable :: integral(:)

    end subroutine TDarkEnergyEqnOfState_SetwTable


    subroutine TDarkEnergyEqnOfState_ReadParams(this, Ini)
    use IniObjects
    use FileUtils
    class(TDarkEnergyEqnOfState) :: this
    class(TIniFile), intent(in) :: Ini
    real(dl), allocatable :: table(:,:)

    end subroutine TDarkEnergyEqnOfState_ReadParams


    subroutine TDarkEnergyEqnOfState_Init(this, State)
    use classes
    class(TDarkEnergyEqnOfState), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    this%is_cosmological_constant = .false.
    this%num_perturb_equations = 0
    
    end subroutine TDarkEnergyEqnOfState_Init

    subroutine TDarkEnergyEqnOfState_PrintFeedback(this, FeedbackLevel)
    class(TDarkEnergyEqnOfState) :: this
    integer, intent(in) :: FeedbackLevel
    end subroutine TDarkEnergyEqnOfState_PrintFeedback
    end module DarkEnergyInterface
