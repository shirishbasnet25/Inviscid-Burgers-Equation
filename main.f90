!Burgers Equation using Different Methods(Inviscid)
program Inviscid_Burger_equation
    implicit none
    Real,dimension(:,:), allocatable::u_n_lax,u_n1_lax,u_n_wen,u_n1_wen,u_up,u_n1_up
    Real,dimension(:), allocatable::x, C_p_up,C_m_up,F_p_up,F_m_up
    Real::F_jp1,F_jm1,F_j,u_jph,u_jmh
    Real::pi = 3.1415927,dx,dt,L=1.0
    Integer::Nt ,Nx = 41,i,j
!Assignment of Value
dx = L / (Nx - 1) ! Step size
Nt = 29; !Number of time steps
!Allocate Variables
Allocate(x(Nx),F_p_up(Nx),F_m_up(Nx),C_p_up(Nx),C_m_up(Nx))
!Generation of spatial points
do i = 1,Nx
    x(i) = (i - 1) * dx
end do

!Allocate variables
Allocate(u_n_lax(Nx,Nt+1), u_n1_lax(Nx,Nt+1),u_n_wen(Nx,Nt+1),u_n1_wen(Nx,Nt+1),u_up(Nx,Nt+1),u_n1_up(Nx,Nt+1))

!Open file for storing the values
open(unit =12,file='output.txt')

!!---------------------------------------------------------------------------------------------
!Inviscid Solution using  Lax Method

!Initial Condition and Boundary Conditions
u_n_lax(:,1) = sin(2 * pi *x(:))
u_n1_lax= u_n_lax
!Writing initial data file
        do i = 1,Nx
            write(12,*)x(i),u_n_lax(i,1);
        end do
        write(12,*)""
        write(12,*)""

    do j = 1,Nt
        dt = dx / maxval(u_n_lax)   !Using stability conditions dt. u_max/ dx <= 1

        do i = 2,Nx-1

        u_n1_lax(i,j) = 0.5 * (u_n_lax(i+1,j) + u_n_lax(i-1,j) ) - 0.25  * dt *  (u_n_lax(i+1,j)**2  -  u_n_lax(i-1,j)**2) / dx;
        end do
        u_n_lax(:,j+1) = u_n1_lax(:,j);

        If (mod(j,5) == 0) then
        !Writing the data into file
        do i = 1,Nx
            write(12,*)x(i),u_n_lax(i,j);
        end do
        write(12,*)""
        write(12,*)""
        end if

    end do
    close(12);
!!-----------------------------------------------------------------------------------------------------------

!!-----------------------------------------------------------------------------------------------------------
!Lax_wendroff method to Invsicid Burgers Equation
open(unit =13, file="output1.txt")
!Initial and Boundary Condition
u_n_wen(:,1) = sin(2 * pi *x(:))
u_n1_wen = u_n_wen
!Writing initial data file
        do i = 1,Nx
            write(13,*)x(i),u_n_wen(i,1);
        end do
        write(13,*)""
        write(13,*)""

do j=1,Nt
    dt = dx / maxval(u_n_wen)
    do i =2,Nx-1
        F_jp1 = 0.5  * u_n_wen(i+1,j)**2
        F_j = 0.5 * u_n_wen(i,j)**2
        F_jm1 =0.5 *  u_n_wen(i-1,j)**2
        u_jph = 0.5 * (u_n_wen(i+1,j) + u_n_wen(i,j))
        u_jmh = 0.5 * (u_n_wen(i,j) + u_n_wen(i-1,j))
        u_n1_wen(i,j) = u_n_wen(i,j) - 0.5 * dt / dx * ( F_jp1 - F_jm1) &
                                + 0.5 * dt**2 / dx**2 * (u_jph * (F_jp1 - F_j) - u_jmh * ( F_j - F_jm1)) ;
    end do
    u_n_wen(:,j+1) = u_n1_wen(:,j);
    If (mod(j,5) == 0) then
        !Writing the data into file
        do i = 1,Nx
            write(13,*)x(i),u_n_wen(i,j);
        end do
        write(13,*)""
        write(13,*)""
        end if
end do
close(13);
!!------------------------------------------------------------------------------------------------
!!--------------------------------------------------------------------------------------------------
!Upwind method for solving Inviscid Burger's Equation
open(unit=14,file="output2.txt")
!Initial and Boundary Conditions
u_up(:,1) = sin(2 * pi *x(:))
u_n1_up = u_up
!Writing initial data file
        do i = 1,Nx
            write(14,*)x(i),u_up(i,1);
        end do
        write(14,*)""
        write(14,*)""
!Determining the value of C_plus,C_minus,F_plus,F_minus
do j=1,Nt
     dt = dx / maxval(u_up)
    do i =1,Nx
        C_p_up(i) = 0.5 * (u_up(i,j) + abs(u_up(i,j)))       !C_plus
        C_m_up(i) = 0.5 * (u_up(i,j) - abs(u_up(i,j)))       !C_minus
        F_p_up(i) =  0.5 * (C_p_up(i) * u_up(i,j))                        !F_plus
        F_m_up(i) = 0.5 * (C_m_up(i) * u_up(i,j))                      !F_minus
    end do

!Applying upwind formula
    do i =2,Nx-1
    u_n1_up (i,j) = u_up(i,j)  - dt / dx * (F_p_up(i) - F_p_up(i-1)) -  dt / dx * (F_m_up(i+1) - F_m_up(i))
    end do
u_up(:,j+1) = u_n1_up(:,j)

 If (mod(j,5) == 0) then
        !Writing the data into file
        do i = 1,Nx
            write(14,*)x(i),u_up(i,j);
        end do
        write(14,*)""
        write(14,*)""
        end if

end do

deallocate(u_n_lax,u_n1_lax,x,u_n1_wen,u_n_wen)
call system('gnuplot -p graph.plt')
!call system('gnuplot -p graph1.plt')
end program




