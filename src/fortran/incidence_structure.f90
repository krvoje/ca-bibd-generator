module incidence_structure
  use utils
  use mtmod

  implicit none

  type IncidenceStructure
     integer, dimension(:,:), allocatable:: incidences
     integer, dimension(:,:), allocatable:: row_intersection
     integer, dimension(:,:), allocatable:: col_intersection
     integer, dimension(:), allocatable:: sum_in_row
     integer, dimension(:), allocatable:: sum_in_col
     integer v
     integer r
     integer k
     integer b
     integer lambda

     ! Heuristic helper vals, computed in the construct method
     integer sum_total
     integer sum_ideal
     integer heuristic_distance
     integer max_heuristic_distance
     integer min_useful_heuristic

     integer generations
  end type IncidenceStructure

contains

  subroutine construct(IS,v,k,lmbd)
    type(IncidenceStructure) IS
    integer v,k,lmbd
    real r_r, b_r
    integer i,j

    is%v=v
    is%k=k
    is%lambda=lmbd

    ! Test for neccessary BIBD conditions
    r_r=is%lambda*(is%v-1)/(is%k-1)
    b_r=r_r*is%v/is%k

    if (r_r/=int(r_r).or.b_r/=int(b_r)) then
       write (*,*) "Invalid BIBD parameters."
       stop
    else
       is%r=int(r_r)
       is%b=int(b_r)
       print *, "(v,k,λ,b,r) = ", is%v,is%k,is%lambda,is%b,is%r
    endif

    allocate(is%incidences(1:is%v,1:is%b))
    allocate(is%row_intersection(1:is%v,1:is%v))
    allocate(is%col_intersection(1:is%b,1:is%b))
    allocate(is%sum_in_row(1:is%v))
    allocate(is%sum_in_col(1:is%b))

    is%incidences(1:is%v,1:is%b)=0
    is%row_intersection(1:is%v,1:is%v)=0
    is%col_intersection(1:is%b,1:is%b)=0
    is%sum_in_row(1:is%v)=0
    is%sum_in_col(1:is%b)=0
    is%sum_total=0

    is%sum_ideal=is%k*is%b
    is%generations=0

    call updateCache(IS)
  end subroutine construct

  subroutine updateCache(IS)
    type(IncidenceStructure) IS

    integer i,j

    do i=1,is%v
       is%sum_in_row(i)=sum(is%incidences(i,:))
    enddo

    do j=1,is%b
       is%sum_in_col(j)=sum(is%incidences(:,j))
    enddo

    do i=1,is%v
       do j=1,is%v
          is%row_intersection(i,j)=dot_product(is%incidences(i,:), is%incidences(j,:))
       enddo
    enddo

    do i=1,is%b
       do j=1,is%b
          is%col_intersection(i,j)=dot_product(is%incidences(:,i), is%incidences(:,j))
       enddo
    enddo

    is%sum_total=sum(is%sum_in_row(1:is%v))

    call calculateHeuristicDistance(IS)
  end subroutine updateCache

  subroutine deconstruct(IS)
    type(IncidenceStructure) IS
    deallocate(is%incidences)
    deallocate(is%row_intersection)
    deallocate(is%col_intersection)
    deallocate(is%sum_in_row)
    deallocate(is%sum_in_col)
  end subroutine deconstruct

  logical function active(IS,row,col)
    type(IncidenceStructure) IS
    integer row,col

    if(is%incidences(row,col) /= 0) then
       active=.True.
    else
       active=.False.
    endif

    return
  end function active

  logical function dormant(IS,row,col)
    type(IncidenceStructure) IS
    integer row,col

    if(is%incidences(row,col) /= 0) then
       dormant=.False.
    else
       dormant=.true.
    endif

    return
  end function dormant

  ! Flips a CA's cell value, and updates the sums and products
  subroutine flip (IS,row,col)

    type(IncidenceStructure) IS
    integer row,col
    integer otherRow,otherCol,newVal

    newVal=abs(is%incidences(row,col)-1)
    is%incidences(row,col)=newVal
    if (newVal==0) then
       call decrement(is%sum_total,1)
       call decrement(is%sum_in_row(row),1)
       call decrement(is%sum_in_col(col),1)
       call decrement(is%row_intersection(row,row),1)
       call decrement(is%col_intersection(col,col),1)
    endif
    if (newVal==1) then
       call increment(is%sum_total,1)
       call increment(is%sum_in_row(row),1)
       call increment(is%sum_in_col(col),1)
       call increment(is%row_intersection(row,row),1)
       call increment(is%col_intersection(col,col),1)
    endif

    do otherRow=1,is%v
       if(is%incidences(otherRow,col)==1) then
          if(newVal==0) then
             call decrement(is%row_intersection(otherRow,row),1)
             call decrement(is%row_intersection(row,otherRow),1)
          endif
          if(newVal==1) then
             call increment(is%row_intersection(otherRow,row),1)
             call increment(is%row_intersection(row,otherRow),1)
          endif
       endif
    enddo

    do otherCol=1,is%b
       if(is%incidences(row,otherCol)==1) then
          if(newVal==0) then
             call decrement(is%col_intersection(otherCol,col),1)
             call decrement(is%col_intersection(col,otherCol),1)
          endif
          if(newVal==1) then
             call increment(is%col_intersection(otherCol,col),1)
             call increment(is%col_intersection(col,otherCol),1)
          endif
       endif
    enddo

    call calculateHeuristicDistance(IS)
  end subroutine flip

  subroutine calculateHeuristicDistance(IS)
    type(IncidenceStructure) IS

    integer i,j

    is%heuristic_distance = 0
    
    do i=1,is%b
       call increment(is%heuristic_distance,abs(is%sum_in_col(i) - is%k))
    enddo
    do i=1,is%v
       do j=1,is%v
          if(i==j) cycle
          call increment(is%heuristic_distance, abs(is%row_intersection(i,j) - is%lambda))
       enddo
    enddo
    is%max_heuristic_distance = is%b*(is%v - is%k) + (is%v**2 - is%v)*(is%r - is%lambda)
  end subroutine calculateHeuristicDistance

  !!***
  !! Checks if the incidence matrix with the other parameters is a BIBD
  !! Returns .True. if this is a 2-(v,k,λ) BIBD
  !! .False. otherwise
  logical function isBIBD(IS)

    type(IncidenceStructure) IS

    integer i,j

    ! if(is%heuristic_distance <= is%max_heuristic_distance ) then
    !      print *, "hd: ", is%max_heuristic_distance, "/", is%heuristic_distance
    !      call writeMatrix(IS)
    ! endif

    if(is%sum_total /= is%sum_ideal) then
       isBIBD=.False.
       !print *, "Sum is non-ideal", is%sum_total, is%sum_ideal
       return
    endif

    ! If the incidence matrix has a row with sum non-equal to r, this is not a BIBD
    do i=1,is%v
       if(is%sum_in_row(i)/=is%r) then
          isBIBD=.False.
          !print *, "Sum in row is non-ideal", is%sum_in_row(i), is%r
          return
       endif
    enddo

    ! If the incidence matrix has a row with sum non-equal to k, this is not a BIBD
    do i=1,is%b
       if(is%sum_in_col(i)/=is%k) then
          isBIBD=.False.
          !print *, "Sum in col is non-ideal", is%sum_in_col(i), is%k
          return
       endif
    enddo

    ! If there are two is%b (rows) whose intersection does not contain lambda is%v,
    ! this is not a BIBD
    do i=1,(is%v-1)
       do j=(i+1),is%v
          if (is%row_intersection(i,j)/=is%lambda) then
             isBIBD=.False.
             !print *, "Intersection is non-lambda", is%row_intersection(i,j)
             return
          endif
       enddo
    end do

    ! Otherwise, we got a BIBD
    isBIBD=.True.
    print *, "Generations: ", is%generations
    return
  end function isBIBD

  subroutine writeMatrix(IS)
    type(IncidenceStructure) IS
    integer row,col

    write(*,*),"["
    do row=1,is%v
       do col=1,is%b
          write (*,"(I1)",advance='no') is%incidences(row,col)
       enddo
       write (*,"(A1)",advance='no') "-"
       write (*,"(I1)",advance='no') is%sum_in_row(row)
       print *
    enddo
    do col=1,is%b
       write (*,"(A1)",advance='no') "|"
    enddo
    print *
    do col=1,is%b
       write (*,"(I1)",advance='no') is%sum_in_col(col)
    enddo
    print *
    write(*,*),"]"
  end subroutine writeMatrix
end module incidence_structure
