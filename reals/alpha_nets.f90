!
!    
! -------------------------------------------------------------------------
!                              alpha_nets.f90
! --------------------------------------------------------------------------
!
! Monte Carlo simulator of a system of spins that can take 3 values 
! corresponding to the vectors (-1,0),(0,alpha) and (1,0) 
! respectively, where alpha is a real parameter. Energy is defined as the
! negative sum of dot products over neighbours with nodes embedded on a network
! taken from a input file.
! From a random initial configuration the program proposes sequencial changes
! and accepts those that minimize the energy following a Boltzmann distribution
! at a given range of temperatures.
! Averages over 'imctot' final states and variances of the energy and the
! proportions of each possible value are calculated, written to a file 
! and printed on the screen for every temp.
 
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
      module topological
      implicit none
      contains
!  DEGREE DISTRIBUTION and CUMULATIVE DEGREE DISTRIBUTION ---
! -----------------------------------------------------------------------------
        subroutine deg_distrib(N, D, max_k, hist_k, cumul_hist_k)
          integer :: i, N, max_k
          integer, dimension(1:N) :: D
          real, dimension(1:max_k) :: hist_k, cumul_hist_k
          real :: nodes
          nodes = real(N)
!          max_k = maxval(D)
          do i = 1, max_k
            hist_k(i) = 0
            cumul_hist_k(i) = 0
          enddo
          do i = 1, N
            hist_k(D(i)) = hist_k(D(i)) + 1
          enddo
          cumul_hist_k(1) = hist_k(1)
          do i = 2, max_k
            cumul_hist_k(i) = cumul_hist_k(i-1) + hist_k(i)
          enddo
          hist_k = hist_k/nodes ! Normalize
          cumul_hist_k = cumul_hist_k/nodes
        end subroutine deg_distrib
!
! AVERAGE NEAREST NEIGBOURS DEGREE
! -----------------------------------------------------------------------------
        subroutine ave_near_deg(N, E, D, V, Pini, Pfin, max_k, k_nn)
          integer :: i, N, E, max_k, k, neigh
          integer, dimension(1:N) :: D, Pini, Pfin
          integer, dimension(1:2*E) :: V
          real, dimension(1:max_k) :: k_nn, hista_k, cumul_hista_k
          real :: nodes, ave_neigh
          nodes = real(N)
!          max_k = maxval(D)
          call deg_distrib(N, D, max_k, hista_k, cumul_hista_k)
          do k = 1, max_k
            k_nn(k) = 0.0
          enddo
          do i = 1,N  ! For each node
            k = D(i)  ! the degree k of node i
            neigh = Pini(i)  ! Starting from the first neighbour of node i
            ave_neigh = 0.0
            do while (neigh < Pfin(i))
              ave_neigh = ave_neigh + D(V(neigh))  ! add the degree of each neighbor of i
              neigh = neigh + 1   ! Go to the following neighbour until it reaches the last one marked by Pfin(i)
            enddo
            ave_neigh = ave_neigh/(float(k)*hista_k(k))
            k_nn(k) = k_nn(k) + ave_neigh
          enddo
!
        end subroutine ave_near_deg 

! CLUSTERING COEFFICIENT
! -----------------------------------------------------------------------------
        subroutine clust_coeff(N, E, D, Pini, Pfin, max_k, clust)
          integer :: i, N, E, max_k, k, T, one, two, check_one 
          integer, dimension(1:N) :: D, Pini, Pfin
          integer, dimension(1:2*E) :: V
          real, dimension(1:max_k) :: c 
          real :: nodes, clust 
          nodes = real(N)
!          max_k = maxval(D)
          do k = 1, max_k
            c(k) = 0.0
          enddo
!
          clust = 0.0
          do i = 1, N
            k = D(i)  ! the degree k of node i      
            T = 0
            do one = Pini(i), Pfin(i)-1  ! Loop over the first neighbour
              do two = one+1, Pfin(i) !  Loop over the second neighbour 
                do check_one = Pini(V(one)), Pfin(V(one)) ! Loop to check if two is in the neigbourhood of one
                  if (V(check_one) == V(two)) then
                    T = T + 1
                  endif
                enddo
              enddo
            enddo
!            write(*,*) T, k
            c(k) = c(k) + (2.0*real(T))/(real(k)*(real(k)-1))
!            write(*,*) c(k)
          enddo
          do k = 1, max_k
            clust = clust + c(k)
          enddo
!          write(*,*) c
          clust = clust/nodes
        end subroutine clust_coeff 
!
      end module topological
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
      module dynamics 
      implicit none
      contains
!
!  FIND INDEX  -----------------------------------------------------------
        subroutine findindex(length, array, val, ini, fin, my_index)
          integer :: my_index, ini, fin, i, val, length
          integer, dimension (1:length) :: array
!
          do i = ini, fin
            if (array(i) == val) then
              my_index = i
              return
            else
              my_index = 0
            endif
          enddo
          return
        end subroutine findindex
!
! SUBROUTINE ENERGY
! -----------------------------------------------------------------------------------
        subroutine energy_function(N,L,alpha, nodes_vec, V, Pini, Pfin,last0,last_1, ene) 
          integer :: summ_1, summ1, i, j, N, L
          real :: summ_bis, alpha
          real(16) :: ene
          integer, dimension(1:N) :: nodes_vec, Pini, Pfin
          integer, dimension(1:L) :: V
          integer :: index_neighbor, last0, last_1
!          integer :: my_index
          ene = 0.0
          summ_1 = 0
          summ_bis = 0.0
          summ1 = 0
          do i = 1, last_1
            do j = Pini(nodes_vec(i)), Pfin(nodes_vec(i))
              call findindex(N, nodes_vec, V(j), 1, N, index_neighbor)
              if (index_neighbor <= last_1) then
                summ_1 = summ_1 - 1  ! criteria equals subtract
              else if (index_neighbor > last0) then
                summ_1 = summ_1 + 1
              endif
            enddo
          enddo ! -------------- contribution of -1 spins
          ene = ene + float(summ_1)
          do i = last_1+1, last0
            do j = Pini(nodes_vec(i)), Pfin(nodes_vec(i))
              call findindex(N, nodes_vec, V(j), 1, N, index_neighbor)
              if (last_1 < index_neighbor .and. index_neighbor <= last0) then
                summ_bis = summ_bis + alpha ! careful, here criteria equals sum*
              endif
            enddo
          enddo
          ene = ene - (alpha*summ_bis)! -------------- contribution of 0 spins (*compensate criteria by substracting)
          do i = last0+1, N
            do j = Pini(nodes_vec(i)), Pfin(nodes_vec(i))
              call findindex(N, nodes_vec, V(j), 1, N, index_neighbor)
              if (index_neighbor <= last_1) then
                summ1 = summ1 - 1 ! careful, here criteria equals sum*
              else if (index_neighbor > last0) then
                summ1 = summ1 + 1
              endif
            enddo
          enddo 
          ene = ene - float(summ1)! -------------- contribution of +1 spins (*compensate criteria by substracting)
! -----------------------------------------------------------------------------
          ene = 0.5*ene  ! ----- count pairs only once
          return
        end subroutine energy_function ! --------------------------------------
!
! -----------------------------------------------------------------------------
!
        subroutine energy_change(N,L, nodes_vec,V, a, Pini, Pfin, init, &
                                 init_bis, fin, fin_bis, diff, &
                                 last0, last_1) 
          integer :: summ, j, N,L,a, index_neighbor, last0, last_1
          real :: summ_bis, e_initial, e_final
          real(16) :: diff
          real :: init, init_bis, fin, fin_bis
          integer, dimension(1:N) :: nodes_vec, Pini, Pfin
          integer, dimension(1:L) :: V
          summ = 0
          summ_bis = 0.0
          do j = Pini(nodes_vec(a)), Pfin(nodes_vec(a))  ! Neighbours of the selected spin
!            write(*,*) 'neigh = ', V(j)
            call findindex(N, nodes_vec,V(j), 1, N,index_neighbor)
            if (index_neighbor <= last_1) then
              summ = summ - 1
            else if (index_neighbor <= last0) then
              summ_bis = summ_bis + alpha
            else
              summ = summ + 1
            endif
!            write(*,*) 'SUMM = ', summ
          enddo 
!          write(*,*) 'a =', a, 'summ = ', summ, 'summ_bis = ',summ_bis
          e_initial = init*summ + init_bis*summ_bis
          e_final = fin*summ + fin_bis*summ_bis
          diff = e_initial - e_final  ! (changing sign to avoid the minus in hamiltonian)
!          write(*,*) 'diff = ', diff
          return
        end subroutine energy_change! --------------------------------------
!!
             
!  MOVE Node Right SUBROUTINE -----------------------------------------------------------
! -----------------------------------------------------------------------------
        subroutine move_node_right(length, array, selected, last)
          integer ::  length, old_last, selected, last, old
          integer, dimension(1:length) :: array 
          old_last = array(last)
          old = array(selected)
          array(last) = old! Switch nodes
          array(selected) = old_last
          last = last - 1
        return
        end subroutine move_node_right
!
!  MOVE Node Left SUBROUTINE -----------------------------------------------------------
! -----------------------------------------------------------------------------
        subroutine move_node_left(length, array, selected, last)
          integer ::  length, old_first, selected, last, old
          integer, dimension(1:length) :: array
          old_first = array(last+1)
          old = array(selected)
          array(last+1) = old
          array(selected) = old_first
          last = last + 1
        return
        end subroutine move_node_left
!
      end module dynamics 
! -----------------------------------------------------------------------------
!! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
!      MAIN
! -----------------------------------------------------------------------------
      program alpha_nets
      use topological
      use dynamics
!      
      implicit none
!
      integer :: spins, ind 
      integer :: last0, last_1 
      integer, dimension(:), allocatable :: V, Pini, Pfin, D
      integer, dimension(:), allocatable :: nodes_vec 
      integer, dimension(:), allocatable :: vec1, vec0, vec_1
      real :: edges, ave_k, nodes,links 
      real :: TIME1,TIME2
      real :: temp, temp0, tempf, temp_step, temp_low, alpha, alpha2, imctot 
      integer :: i,j,t,k,N,L,E,istep,iseed
      integer :: MCtot,steps_low, steps, imc, temp_tot, reps, reps_tot
      real(16) :: DE,energyt,energy2,energy
      integer :: num0, num1, num_1, numx, numx_abs
      integer :: new_spin, new_num0, new_num1, new_num_1
      integer(16) :: num0total, num1total, num_1total, numxtotal, numx_abstotal
      integer(16) :: N2, num02, num12, num_12, numx2 
      integer :: a
      integer(2) :: b
      real(16) :: f_energy, f_num0, f_num1, f_num_1, f_numx, f_numx_abs
      real(16) :: f_energy2, f_num02, f_num12, f_num_12, f_numx2  ! Squares
      real(16) :: vare, var0, var1, var_1, varx  ! Variances 
      real :: r1279, randunif
      character(15):: output
      character(12) :: input
      character(30):: date
      logical :: left_zero, left_right, zero_left, zero_right, right_zero
      real :: initial_spin, final_spin, initial_spin_bis, final_spin_bis
!      integer :: my_index 
!
      NAMELIST/DADES/output,input,iseed,alpha,temp0,tempf,temp_tot,&
                     temp_low,steps_low, steps, reps_tot ! Initial data
      OPEN(UNIT=10,FILE='alpha3states_nets.dat')        ! Input
      READ(10,DADES)
      CLOSE(10)
!
      OPEN(UNIT=13,FILE=output//".csv")  ! Results file
!
! ----------------------------------------------------------------------!
!      Read the file and store N and E
! -----------------------------------------------------------------------------
!
      N = 0
      E = 0
      open(unit=15,file=input//".net")
      do
        read(15,*,end=20) a, b
        E = E + 1
        if (a > N) then
          N = a
        endif
        if (b > N) then
          N = b
        endif
      enddo
!
  20  close(15)
      edges = float(E)
      nodes = float(N)
      links = 2.0*edges
      ave_k = links/nodes
      L = int(links)
      write(*,*) 'Number of nodes = ', N
      write(*,*) 'Number of edges = ', E
      write(*,*) 'Average degree <k> = ', ave_k
      allocate(V(1:L))  ! Neigbors vector
      allocate(Pfin(1:N))  ! final neigbour pointer to V
      allocate(Pini(1:N))  ! initial neigbour pointer to V
      allocate(D(1:N))    ! Degrees vector
      allocate(nodes_vec(1:N))
      allocate(vec1(1:N))
      allocate(vec0(1:N))
      allocate(vec_1(1:N))
!
      temp0 = temp0*ave_k  ! Rescaling noise by J (constant coupling)
      tempf = tempf*ave_k
      temp_low = temp_low*ave_k
!
!  INITIALIZE VECTORS AND VARIABLES -----------------------------------------
      do i = 1, N
        D(i) = 0
        Pini = 0
        Pfin = 0
      enddo
      do i = 1, L
        V(i) = 0
      enddo
      do i = 1, N
        nodes_vec(i) = i
      enddo
!
!  READ THE FILE, STORE DEGREES INFORMATION ---------------------------------
      open(unit=15,file=input//".net")
      do i = 1, E
        read(15,*) a, b
        ind = a
        D(ind) = D(ind) + 1
        ind = b
        D(ind) = D(ind) + 1
      enddo
      close(15)
!      write(*,*) D
!
!  INITIALIZE Pini AND Pfin  ------------------------------------------------
      Pini(1) = 1
      do i = 2, N
        Pini(i) = Pini(i-1) + D(i-1)
        Pfin(i-1) = Pini(i) - 1  ! Pfin[i-1] and Pini[i] must be consecutive
      enddo
      Pfin(N) = 2*E
!
!  READ THE FILE AND STORE NEIGHBOURS MATRIX IN V UPDATE  -------------------
      j = 1
      do i = 1, N
        open(unit=15,file=input//".net")
        do k = 1, E
          read(15,*) a, b
          ind = a
          if (ind == i) then
            V(j) = b
            j = j + 1
          endif
          ind = b
          if (ind == i) then
            V(j) = a
            j = j + 1
          endif
        enddo
        close(15)
      enddo
!
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
      write(*,*) 'Number of spins = ', N
!      write(*,*) 'MCS = ', MCtot
      temp_step = (tempf - temp0)/float(temp_tot)
      N2 = N*N
      alpha2 = alpha*alpha
!
! -------------------------------------------------------------------
      CALL setr1279(iseed)    
      temp = temp0
!
      CALL CPU_TIME(TIME1)          ! Time control
      do t = 1, temp_tot !  TEMPERATURE LOOP ---------------------------------------------
        if (temp < temp_low) then
          MCtot = steps_low
        else
          MCtot = steps
        endif
        energyt=0
        num0total=0
        num1total=0
        num_1total=0
        numxtotal=0
        numx_abstotal=0
        energy2=0
        num02=0
        num12=0
        num_12=0
        numx2=0
        do reps = 1, reps_tot  !     REPETITIONS LOOP
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------
! --- INITIAL CONFIGURATION. -----------------------------------------------
!  Random spins -1, 0 or +1 ----------------------- 
!
!  All spins are supposed to start in -1, if the random generator sets a 0
! the spin will be moved like in a SIR model to the medium group and if 
! it sets the spin to +1 the node will be moved twice in a row to the right 
! group
! --------------------------------------------------------------------------- 
!
          num_1 = 0; num0 = 0; num1 = 0  ! Number of agents in each state
          do i = 1, N
            spins = mod(int(3*r1279()),3) - 1     
!            write(*,*) spins 
            if (spins == -1) then
              num_1 = num_1 + 1
              vec_1(num_1) = i  ! Vector that contains node labels of nodes in state -1
            else if (spins == 0) then
              num0 = num0 + 1
              vec0(num0) = i    ! Vector that contains node labels of nodes in state 0
            else
              num1 = num1 + 1 
              vec1(num1) = i     ! Vector that contains node labels of nodes in state +1
            endif
          end do  ! --------------------------------------------------
          do i = 1, num_1
            nodes_vec(i) = vec_1(i)  ! Vector that contains node labels of all nodes, ordered in 3 boxes (-1, 0 and +1), separated by the pointers last_1 and last0
          enddo
          do i = 1, num0
            nodes_vec(num_1 + i) = vec0(i)
          enddo
          do i = 1, num1
            nodes_vec(num_1+num0 + i) = vec1(i)
          enddo
!          write(*,*) 'nodes_vec =', nodes_vec
!
!   INITIALIZE Pointers last_1 and last0
!    
          last_1 =  num_1  ! Points to the last node in state -1
          last0 = num_1 + num0      ! Points to the last node in state 0
!          write(*,*) 'last_1 = ', last_1, 'last0 =' , last0
! ----------------------------------------------------------------------------------------------------------
!   INITIAL ENERGY
!          
          call energy_function(N,L,alpha, nodes_vec, V, Pini, Pfin, last0, last_1,energy)          
!          write(*,*) 'Initial energy =', energy
          new_num_1 = num_1
          new_num0 = num0
          new_num1 = num1
! -------------------------------------------------------------------
!     Metropolis. GLAUBER Dynamics
! -------------------------------------------------------------------
! Markov chain -------------------------------------------------------------------------------------------           
          DO imc = 1,MCtot
            DO ISTEP = 1,N         ! Each MC has N flip attempts
              left_zero = .false.; left_right=.false.; zero_left=.false. 
              zero_right=.false.; right_zero=.false.
              a = mod(int(N*r1279()), N) + 1  ! Pick a random spin
              b = mod(int(2*r1279()), 2)   ! Change proposal
!              call findindex(N, nodes_vec, a, 1, N, index_a)
              if (a <= last_1) then
                initial_spin = -1
                initial_spin_bis = 0  ! Second component of the initial vector -----------------------
                if (b == 0) then
                  final_spin = 0.0 
                  final_spin_bis = alpha
                  left_zero = .true.                
!                  write(*,*) 'left_zero ----------------------'
                  call energy_change(N,L, nodes_vec,V, a, Pini, Pfin, &
                                  initial_spin, initial_spin_bis,&
                                  final_spin,final_spin_bis, DE,&
                                  last0, last_1)
                  new_spin = 0              !  -1 to 0
                  new_num1=num1
                  new_num0=num0+1
                  new_num_1=num_1-1
                else
                  final_spin = 1.0
                  final_spin_bis = 0.0
                  left_right = .true.
!                  write(*,*) 'left_right ---------------------'
                  call energy_change(N,L, nodes_vec,V,a, Pini, Pfin, &
                                  initial_spin, initial_spin_bis,& 
                                  final_spin,final_spin_bis, DE,&
                                  last0, last_1)
                  new_spin = 1              !  -1 to +1
                  new_num1=num1+1
                  new_num0=num0
                  new_num_1=num_1-1
                endif
              else if (a <= last0) then
                initial_spin = 0
                initial_spin_bis = alpha  ! Second component of the initial vector ------------------
                if (b == 0) then
                  final_spin = -1.0
                  final_spin_bis = 0.0
                  zero_left = .true.
!                  write(*,*) 'zero_left ---------------------'
                  call energy_change(N,L, nodes_vec,V, a, Pini, Pfin, &
                                  initial_spin, initial_spin_bis,&
                                  final_spin,final_spin_bis, DE,&
                                  last0, last_1)
                  new_spin = -1             ! 0 to -1
                  new_num1=num1
                  new_num0=num0-1
                  new_num_1=num_1+1
                else
                  final_spin = 1.0
                  final_spin_bis = 0.0
                  zero_right = .true.
!                  write(*,*) 'zero_right ---------------------'
                  call energy_change(N,L, nodes_vec,V, a,Pini, Pfin, &
                                  initial_spin, initial_spin_bis,&
                                  final_spin,final_spin_bis, DE,&
                                  last0, last_1)

                  new_spin = 1              !  0 to +1
                  new_num0=num0-1
                  new_num1=num1+1
                  new_num_1=num_1
                endif
              else
                initial_spin = 1.0
                initial_spin_bis = 0.0  ! Second component of the opinion vector
                if (b == 0) then
                  final_spin = 0.0
                  final_spin_bis = alpha
                  right_zero = .true.
!                  write(*,*) 'right_zero ---------------------'
                  call energy_change(N,L, nodes_vec,V, a , Pini, Pfin, &
                                  initial_spin, initial_spin_bis,&
                                  final_spin,final_spin_bis, DE,&
                                  last0, last_1)
                  new_spin = 0              ! +1 to 0
                  new_num0=num0+1
                  new_num1=num1-1
                  new_num_1=num_1
                else
                  final_spin = -1.0
                  final_spin_bis = 0.0
!                  write(*,*) 'right_left ---------------------'
                  call energy_change(N,L, nodes_vec,V, a, Pini, Pfin, &
                                  initial_spin, initial_spin_bis,&
                                  final_spin,final_spin_bis, DE,&
                                  last0, last_1)
                  new_spin = -1             ! +1 to -1
                  new_num1=num1-1
                  new_num_1=num_1+1
                  new_num0=num0
                endif
              endif
!      Accept or reject ------------------------------------------------------------------------------
              if (DE < 0.0) then ! ACCEPT ---------------------------------------------------------
                energy = energy + DE ! -------------------------------------------------------------
                if (left_zero .eqv. .true.) then  ! -1 to 0
                  call move_node_right(N, nodes_vec, a, last_1)
                else if (left_right .eqv. .true.) then ! -1 to +1
                  call move_node_right(N, nodes_vec, a, last_1) ! passes through 0
                  call move_node_right(N, nodes_vec, last_1+1, last0)
                else if (zero_left .eqv. .true.) then  ! 0 to -1 
                  call move_node_left(N, nodes_vec, a, last_1)
                else if (zero_right .eqv. .true.) then  ! 0 to +1
                  call move_node_right(N, nodes_vec, a, last0)
                else if (right_zero .eqv. .true.) then ! +1 to 0
                  call move_node_left(N, nodes_vec, a, last0)
                else
                  call move_node_left(N, nodes_vec, a, last0)  ! passes through 0
                  call move_node_left(N, nodes_vec, last0, last_1)
                endif
                num0=new_num0
                num1=new_num1
                num_1=new_num_1
!                write(*,*) 'a = ', a, 'ACCEPT'
!                write(*,*) 'nodes_vec = ', nodes_vec
!                write(*,*) 'last_1 = ', last_1, 'last0 =' , last0
!                write(*,*) 'energy = ', energy
              else
                randunif = r1279() 
                if (randunif < exp(-DE/temp)) then   ! ACCEPT ------------------------------
                  energy = energy + DE ! --------------------------------------------------------
                if (left_zero .eqv. .true.) then  ! -1 to 0
                  call move_node_right(N, nodes_vec, a, last_1)
                else if (left_right .eqv. .true.) then ! -1 to +1
                  call move_node_right(N, nodes_vec, a, last_1) ! passes through 0
                  call move_node_right(N, nodes_vec, last_1+1, last0)
                else if (zero_left .eqv. .true.) then  ! 0 to -1 
                  call move_node_left(N, nodes_vec, a, last_1)
                else if (zero_right .eqv. .true.) then  ! 0 to +1
                  call move_node_right(N, nodes_vec, a, last0)
                else if (right_zero .eqv. .true.) then ! +1 to 0
                  call move_node_left(N, nodes_vec, a, last0)
                else
                  call move_node_left(N, nodes_vec, a, last0)  ! passes through 0
                  call move_node_left(N, nodes_vec, last0, last_1)
                endif
                  num0=new_num0
                  num1=new_num1
                  num_1=new_num_1
!                  write(*,*) 'a = ', a, 'ACCEPT'
!                  write(*,*) 'nodes_vec = ', nodes_vec
!                  write(*,*) 'last_1 = ', last_1, 'last0 =' , last0
!                  write(*,*) 'energy = ', energy
                endif
              endif
!
            end do  ! Ends steps (1 to N) loop ------------------------------
          end do  ! ends imc loop
! 
          energyt=energyt+energy             ! Summation
          num0total=num0total+num0
          num1total=num1total+num1
          num_1total=num_1total+num_1
          numx=num1-num_1
          numxtotal=numxtotal+numx
          numx_abstotal=numx_abstotal+abs(numx)
          energy2=energy2+energy*energy
          num02=num02+num0*num0
          num12=num12+num1*num1
          num_12=num_12+num_1*num_1
          numx2=numx2+numx*numx
        end do   ! ends reps loop
!     
        imctot = float(reps_tot)
        f_energy=energyt/imctot              ! Normalization and errors
        f_num0=float(num0total)/imctot
        f_num1=float(num1total)/imctot
        f_num_1=float(num_1total)/imctot
        f_numx=float(numxtotal)/imctot
        f_numx_abs=float(numx_abstotal)/imctot
        f_energy2=energy2/imctot
        f_num02=float(num02)/imctot
        f_num12=float(num12)/imctot
        f_num_12=float(num_12)/imctot
        f_numx2=float(numx2)/imctot
        vare=energy2-(energy*energy)
        var0=num02-(num0*num0)
        var1=num12-(num1*num1)
        var_1=num_12-(num_1*num_1)
        varx=numx2-(numx_abs*numx_abs)
! -------------------------------------------------------------------------
!     WRITING
! ---------------------------------------------------------------------
        write(*,*) N, N2, temp, energy, f_num_1, f_num0, f_num1
        write(13,*) N, N2, temp, f_energy, vare, f_num0, f_num1, f_num_1, &
                    f_numx_abs, var0, var1, var_1, varx
        temp = temp + temp_step
      END DO  ! Ends temp loop
!
      CALL CPU_TIME(TIME2)            ! Time control
      CALL FDATE(DATE)
!      flips = N*MCtot
!      VEL = float(flips)/(TIME2-TIME1)
      WRITE(*,*) DATE
      WRITE(*,*) 'CPUTIME=',TIME2-TIME1
!      WRITE(*,*) 'VELOCITY=', VEL, 'ATTEMPTED FLIPS PER SECOND'
!
      CLOSE(13)
      END program alpha_nets

