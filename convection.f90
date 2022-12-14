program LowPr
    use iso_fortran_env
    use all
    
    implicit none
    include "mpif.h"

    type(everything) :: all
    real(kind=real64) :: time
    integer :: step

    call all%read_inputs()
    call all%initialize()

    time = 0.; step = 0

    do while (time < all%total_time)
        call all%VelocitySolve()
        call all%Timestep()
        time = time + all%dt
        step = step + 1
        if (mod(step, 100) == 1 .or. time>= all%total_time) then
            print*, 'Time, step, maxV : ', time, step, maxval(abs(all%Vmax))
            call all%WriteTField() 
        end if 
    end do

    call all%WriteInfos()

end program LowPr