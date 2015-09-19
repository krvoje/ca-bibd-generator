module tabu_list
    use utils

    implicit none

    type TabuList
        integer, dimension(:), allocatable:: tabus
        integer head, tail, length
    end type TabuList

    contains
    
    subroutine mkTabuList(tl, length)
        type(TabuList) tl
        integer length

        tl%head = 1
        tl%length = length
        allocate(tl%tabus(1:length))
        tl%tabus(1:length)=-1
    end subroutine

    subroutine push(tl, val)
        type(TabuList) tl
        integer length, val

        tl%tabus(tl%head) = val
        call increment(tl%head, 1)
        if(tl%head > tl%length) tl%head = 1
    end subroutine

    subroutine delTabuList(tl)
        type(TabuList) tl
        deallocate(tl%tabus)
    end subroutine

    logical function is_in_list(tl, val) 
        type(TabuList) tl
        integer i, val
        do i = 1, tl%length
            if(tl%tabus(i) .eq. val) then 
                is_in_list = .True.
                return
            endif
        enddo        
        is_in_list=.False.
        return
    end function

end module tabu_list
