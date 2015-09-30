module n_opt
  implicit none
  public::nOpt
  
contains
  
  recursive logical function nOpt(n,topOpt,is) result (successfulOpt)
    use incidence_structure
    implicit none
    type(IncidenceStructure) is
    integer n,row,col1,col2,topOpt
    integer col1_delta, col2_delta

    logical sumTotal_low, sumTotal_high
    logical sumInRow_low, sumInRow_high

    successfulOpt=.false.
    
    ! Skip criteria for the n-opt:
    ! We reached the end of the n-opt descent, return
    if(n<1) return
    if(abs(is%heuristic_distance) > is%max_heuristic_distance*n / 2) return
    ! For each row try exchanging each two opposite-state cells
    do row=1,is%v
       do col1=1,is%b - 1
          do col2=col1 + 1,is%b
             ! Do not exchange cells with the same state
             if(dormant(is,row,col1).and.dormant(is,row,col2)) cycle
             if(active(is,row,col1).and.active(is,row,col2)) cycle

             if(abs(is%sum_in_col(col1)-is%k) > n) cycle
             if(abs(is%sum_in_col(col2)-is%k) > n) cycle
             
             ! Exchange the cells' state             
             call flip(is,row,col1)
             call flip(is,row,col2)
          
             ! Try testing for BIBD
             if(isBIBD(is)) then
                print *, topOpt, "-opt yielded a BIBD"
                successfulOpt=.true.
                return
             ! Else try descending
             else if(nOpt(n-1,topOpt,is)) then
                successfulOpt=.true.
                return
                ! If the descent failed, flip the cell back and continue
             else
                call flip(is,row,col1)
                call flip(is,row,col2)
             endif
          enddo
       enddo
    enddo
    ! Upon reaching here we return .FALSE.
    return
  end function nOpt

end module n_opt
