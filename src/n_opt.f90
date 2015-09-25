module n_opt
  implicit none
  public::nOpt
  
contains
  
  recursive logical function nOpt(n,topOpt,is) result (successfulOpt)
    use incidence_structure
    implicit none
    type(IncidenceStructure) is
    integer n,row,col,topOpt  

    logical sumTotal_low, sumTotal_high
    logical sumInRow_low, sumInRow_high

    successfulOpt=.false.
    
    ! Skip criteria for the n-opt:
    ! We reached the end of the n-opt descent, return
    if(n<1) return
    ! If it does not make sence to attempt the optimization, return
    if (abs(is%heuristic_distance) > is%max_heuristic_distance * n) return 
    ! If it does not make sence to attempt the optimization for this
    ! number of steps, try rerunning the opt algorithm for a lower step.
    if (abs(is%sum_ideal - is%sum_total) /= n) then
       successfulOpt = nOpt(n-1,n-1,is)
       return
    endif

    ! Trying to speed things up by calculating the flags beforehand
    sumTotal_low = (is%sum_total <= is%sum_ideal - n)
    sumTotal_high = (is%sum_total >= is%sum_ideal + n)

    ! For each cell in the incidence matrix try optimizing
    do row=1,is%v
       ! Trying to speed things up by calculating the flags beforehand
       sumInRow_low = (is%sum_in_row(row) <= is%r - n)
       sumInRow_high = (is%sum_in_row(row) >= is%r + n)
       do col=1,is%b
          ! No point in flipping a cell state, the total sum is too
          ! far away from the ideal sum
          if(dormant(is,row,col)) then
             if(.not.sumTotal_low) cycle
             if(.not.sumInRow_low .and. is%sum_in_col(col) > is%k - n) cycle
          endif
          if(active(is,row,col)) then
             if(.not.sumTotal_high) cycle
             if(.not.(sumInRow_high) .and. is%sum_in_col(col) < is%k + n) cycle
          endif

          ! Flip the cell, and update the flags:
          call flip(is,row,col)
          sumInRow_low = (is%sum_in_row(row) <= is%r-n)
          sumInRow_high = (is%sum_in_row(row) >= is%r+n)
          
          ! If this is our final step, then try testing for BIBD          
          if(n == 1) then
             if (isBIBD(is)) then
                print *, topOpt, "-opt yielded a BIBD"
                successfulOpt=.true.
                return
             endif
          
          ! Otherwise try descending
          else if(nOpt(n-1,topOpt,is)) then
             successfulOpt=.true.
             return

          ! If the descent failed, flip the cell back and fail
          else
             call flip(is,row,col)
             successfulOpt=.false.
             return
          endif
       enddo
    enddo
    ! Upon reaching here we return .FALSE.
    return
  end function nOpt

end module n_opt
