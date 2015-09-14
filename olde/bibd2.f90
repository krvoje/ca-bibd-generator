
    program bibd

        implicit none

        logical isBIBD
        real nasumice

        integer, dimension(:,:), allocatable:: A
        integer v,k,b,lmbd,r,i,j,stay_awake,stay_asleep
        real r_r,b_r
        character(100) f

        if (command_argument_count() /= 3) then
            write (*,*) "Upotreba:"
            call get_command_argument(0,f)
            write (*,*) f ,"v k λ tabulen"
            stop
        else
            call get_command_argument(1,f)
            read (f,"(I10)") v
            call get_command_argument(2,f)
            read (f,"(I10)") k
            call get_command_argument(3,f)
            read (f,"(I10)") lmbd
        endif

        r_r=lmbd*(v-1)/(k-1)
        b_r=r_r*v/k

        stay_awake=0
        stay_asleep=0

        if (r_r/=int(r_r).or.b_r/=int(b_r)) then
            write (*,*) "Ne postoji BIBD s tim parametrima."
            stop
        else
            r=int(r_r)
            b=int(b_r)

            print *, "v,k,λ,b,r=", v,k,lmbd,b,r
        endif

        call posij()
        allocate(A(1:v,1:b))

        a(1:v,1:b)=int(nasumice())

        do i=1,v
            do j=1,b
                if (a(i,j)/=1 .and. a(i,j)/=0) then
                    print *, "greška u rendom đeneratoreu"
                    print *,a(i,j)
                endif
            enddo
        enddo

        do while(isBIBD(A,v,k,lmbd,b,r).eqv..false.)
            call randomCA_BIBD(A,v,k,lmbd,b,r,stay_awake,stay_asleep)
            print *, "Našli smo jednu matricu incidencije BIBD-a s tim parametrima:"
            do i=1,v
            do j=1,b
                write (*,"(I1)",advance='no') A(i,j)
                enddo
                print *
            enddo
            stop
        enddo

        deallocate(A)
    end program

    logical function isBIBD(A,v,k,lmbd,b,r)
        ! Funkcija vraća .True. ako je A BIBD 2-(v,k,λ)
        ! Inače vraća .False.

        integer r,b,v,k,lmbd,i,j
        integer A(1:v,1:b)

        do i=1,v
            if (sum(A(i,:))/=r) then
                isBIBD=.False.
                return
            endif
        enddo

        do i=1,b
            if (sum(A(:,i))/=k) then
                isBIBD=.False.
                return
            endif
        enddo

        do i=1,(v-1)
            do j=(i+1),v
                if (dot_product(A(i,:),A(j,:))/=lmbd) then
                    isBIBD=.False.
                    return
                endif
            enddo
        end do
        isBIBD=.True.
        write (*,*) "Je BIBD!"
        return
    end function

    subroutine randomCA_BIBD(A,v,k,lmbd,b,r,stay_awake,stay_asleep)
        ! TODO: prepisat is Sagea šta funkcija radi
        use mtmod

        implicit none

        integer r,b,v,k,lmbd,i,j,u_koliko, tocka,suma,faktor,suma_j,suma_i ! suma je suma svih incidencija, računa se u hodu da se uštedi na vremenu i prostoru. faktor je faktor hoće li stanica oživjeti ili ne
        integer A(1:v,1:b),stay_awake,stay_asleep,opt_i,opt_j
        logical Kraj,isBIBD,oneOpt, twoOpt
        integer tabulen
        common tabulen
        real nasumice

        integer optc, lasti,lastj

        kraj=.false.
        suma=0
        faktor=0

        optc=0

        do i=1,v
            suma = suma + sum(A(i,:))
        enddo

        i=int(nasumice()*(v-1))+1
        j=int(nasumice()*(b-1))+1
        lasti=i
        lastj=j

        do while(Kraj.eqv..false.)
            ! TODO: vidjeti, ovdje je bio 2-opt, ali zanemarit jer vjerojatno puno previše napumpava složenost
            if (abs(suma-v*r)==1) then  !ako opće ima smisla raditi 1-opt korak
                optc=optc+1
                if (oneOpt(A,v,k,lmbd,b,r,suma).eqv..true.) then
                    Kraj=.true.
                            print *, "broj 1-opt koraka sve skupa:", optc
                    return
                endif
            endif
            ! Inicijalizacija
            faktor=0

            i=int(nasumice()*(v))+1
            j=int(nasumice()*(b))+1


            lasti=i
            lastj=j


            suma_j=sum(A(:,j))
            suma_i=sum(A(i,:))

            ! Ako je stanica mrtva
            if (A(i,j)==0) then
                !print *,':('
                do tocka=1,v
                    if (tocka==i) cycle
                    if (dot_product(A(i,:),A(tocka,:))<lmbd) faktor=faktor+1 ! Najviše v-1
                enddo
                if (suma<v*r) faktor=faktor+1 ! Najviše 1
                if (suma_j<k) faktor=faktor+1
                if (suma_i<r) faktor=faktor+1
                if(int(nasumice()*(v+2))<faktor) then
                    A(i,j)=1
                    suma=suma+1
                endif
            ! Ako je stanica živa
            else if (A(i,j)==1) then
                !print *,':)'
                do tocka=1,v
                    if (tocka==i) cycle
                    if (dot_product(A(i,:),A(tocka,:))>lmbd) faktor=faktor+1
                enddo
                if (suma>v*r) faktor=faktor+1
                if (suma_j>k) faktor=faktor+1
                if (suma_i>r) faktor=faktor+1
                if(int(nasumice()*(v+2))<faktor) then ! Eksperimenti, dat prednost živim incidencijama
                    A(i,j)=0
                    suma=suma-1
                endif
            endif
        enddo
    end subroutine randomCA_BIBD

    ! Funkcije za random računanje
    subroutine posij()
        use mtmod

        external ZBQLINI

        real nasumice

        call ZBQLINI(time())
        call sgrnd(time())
    end subroutine

    real function nasumice()
        use mtmod

        double precision ZBQLU01,ZBQLUAB, ZBQLEXP
        integer i

        !nasumice=ZBQLU01(i)
        nasumice=grnd() !<- Ova implementacija Mersenne twistera je koma za ovu svrhu
        if (nasumice==1) nasumice=0 ! Računamo da se nikad ne dobije 1, jer bi to potrgalo sve
        return

    end function

    ! Rutine za 1-opt i 2-opt, da bude sve skupa čitljivije:
    logical function oneOpt(A,v,k,lmbd,b,r,suma)
        implicit none

        integer v,k,lmbd,b,r,suma,i,j
        integer A(1:v,1:b)

        logical isBIBD

        oneOpt=.false.
        !print *, "Krećemo u 1-opt!"
        do i=1,v
            do j=1,b
                if(&
                (A(i,j)==0 .and. suma==v*r-1 .and. sum(A(i,:))==r-1 .and. sum(A(:,j))==k-1 )&
                .or.&
                (A(i,j)==1 .and. suma==v*r+1 .and. sum(A(i,:))==r+1 .and. sum(A(:,j))==k+1) ) then
                    A(i,j)=abs(A(i,j)-1)
                    if (isBIBD(A,v,k,lmbd,b,r).eqv..True.) then
                        write (*,*) "Ovaj BIBD smo dobili s 1-opt"
                        oneOpt=.true.
                        return
                    else
                        !write (*,*) "s 1-opt nismo dobili BIBD"
                        A(i,j)=abs(A(i,j)-1)
                        oneOpt=.false.
                        return
                    endif
                endif
            enddo
        enddo
    end function oneOpt


    logical function twoOpt(A,v,k,lmbd,b,r,suma)
        implicit none

        integer v,k,lmbd,b,r,opt_i,opt_j,suma,i,j
        integer A(1:v,1:b)
        logical isBIBD

        twoOpt=.false.
        print *, "Krećemo u 2-opt!"
        do i=1,v
            do j=1,b
                if(&
                (A(i,j)==0 .and. suma==v*r-2 .and. sum(A(i,:))==r-2 .and. sum(A(:,j))==k-2 )&
                .or.&
                (A(i,j)==1 .and. suma==v*r+2 .and. sum(A(i,:))==r+2 .and. sum(A(:,j))==k+2) ) then
                    A(i,j)=abs(A(i,j)-1)
                    do opt_i=1,v ! ! 1-opt podkorak
                        if (opt_i==i) cycle
                        do opt_j=1,v
                            if (opt_j==j) cycle
                            if(&
                            (A(opt_i,opt_j)==0 .and. suma==v*r-1 .and. sum(A(opt_i,:))==r-1 .and. sum(A(:,opt_j))==k-1 )&
                            .or.&
                            (A(opt_i,opt_j)==1 .and. suma==v*r+1 .and. sum(A(opt_i,:))==r+1 .and. sum(A(:,opt_j))==k+1) )&
                            then
                                A(opt_i,opt_j)=abs(A(opt_i,opt_j)-1)
                                if (isBIBD(A,v,k,lmbd,b,r).eqv..True.) then
                                    write (*,*) "Ovaj BIBD smo dobili s 2-opt"
                                    twoOpt=.true.
                                    return
                                else
                                    !write (*,*) "s 2-opt nismo dobili BIBD"
                                    A(opt_i,opt_j)=abs(A(opt_i,opt_j)-1)
                                endif
                            endif
                        enddo
                    enddo
                    ! Ako ono gore nije uspjelo, onda nije ni ovo
                    A(i,j)=abs(A(i,j)-1)
                    twoOpt=.false.
                    return
                endif
            enddo
        enddo
    end function
