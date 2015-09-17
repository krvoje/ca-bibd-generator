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

    if(n<1) return
    if (abs(is%heuristic_distance) > is%max_heuristic_distance * n) return
    if (abs(is%sum_ideal - is%sum_total) /= n) then
       successfulOpt = nOpt(n-1,n-1,is)
       return
    endif

    sumTotal_low = (is%sum_total <= is%sum_ideal - n)
    sumTotal_high = (is%sum_total >= is%sum_ideal + n)

    do row=1,is%v
       sumInRow_low = (is%sum_in_row(row) <= is%r - n)
       sumInRow_high = (is%sum_in_row(row) >= is%r + n)
       do col=1,is%b
          if(dormant(is,row,col)) then
             if(.not.sumTotal_low) cycle
             if(.not.sumInRow_low .and. is%sum_in_col(col) > is%k - n) cycle
          endif
          if(active(is,row,col)) then
             if(.not.sumTotal_high) cycle
             if(.not.(sumInRow_high) .and. is%sum_in_col(col) < is%k + n) cycle
          endif
          call flip(is,row,col)
          sumInRow_low = (is%sum_in_row(row) <= is%r-n)
          sumInRow_high = (is%sum_in_row(row) >= is%r+n)
          if(n == 1) then
             if (isBIBD(is)) then
                print *, topOpt, "-opt yielded a BIBD"
                successfulOpt=.true.
                return
             endif
          else if(nOpt(n-1,topOpt,is)) then
             successfulOpt=.true.
             return
          else
             call flip(is,row,col)
             sumInRow_low = (is%sum_in_row(row) <= is%r-n)
             sumInRow_high = (is%sum_in_row(row) >= is%r+n)
             successfulOpt=.false.
             return
          endif
       enddo
    enddo
    return
  end function nOpt

end module n_opt
