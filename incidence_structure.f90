module incidence_structure
  type IncidenceStructure
     integer, dimension(:,:), allocatable:: incidences
     integer, dimension(:), allocatable:: sumInRow
     integer, dimension(:), allocatable:: sumInCol
     integer v,k,b,lmbd,r,sumTotal
  end type IncidenceStructure
  
end module incidence_structure
