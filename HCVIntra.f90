PROGRAM New_model

  IMPLICIT NONE 
  
  INTEGER, DIMENSION(8) :: rna,syn
  ! rna(1)=cyt 2=ribo 3=ip 4=ids 5=vms 6=p5B 7=ds 8=pack 
  ! syn(1)=cyt 2=ribo 3=ip 4=ids 5=vms 6=p5B 7=ds 8=pack 
  INTEGER, DIMENSION(2) :: ns5b ! 1=cyt 2=vms
  
  INTEGER :: ribo,tc,pp,s,ns34a,ns5a,hf,hfc

  INTEGER :: i,j,k,iseed,moi_v,moi_s
  INTEGER :: nrun,iseed0,synentry
  
  DOUBLE PRECISION :: k_tc,k_trans,k_cleavage ! -- translation rates 
  DOUBLE PRECISION :: k_rp_5b,k_hfc,k_rids,k_repl,k_init  ! -- replication rates 
  DOUBLE PRECISION :: k_rip,k_assembly,k_out_rp,k_s ! -- replication rates 
  DOUBLE PRECISION :: k_sp_5b,k_sids,k_srepl,k_sinit,k_sip,k_out_sp ! -- replication rates for syn
  DOUBLE PRECISION :: k_deg_r,k_deg_syn ! -- degradation rates
  DOUBLE PRECISION, DIMENSION(2) :: k_deg_ns5b 
  DOUBLE PRECISION :: k_deg_si,k_deg_s,k_deg_vms,k_deg_ns,syntime

  DOUBLE PRECISION, DIMENSION(33) :: b 
  DOUBLE PRECISION :: random,btot,time,time_end,r,tau,tout
  DOUBLE PRECISION :: ba,br
        
  iseed = 571723
  
  DO nrun=1,1
  
  iseed0=iseed


  !== readin == !

  moi_v = 1
  moi_s = 0
  synentry = 0
  syntime = 0.0d0
  ba = 1.0d0
  br = 1.0d0

  !== translation rates == ! 

  k_tc =1.0d0
  k_trans = 180.0d0
  k_cleavage = 9.0d0

  !== replication rates == !

  k_rp_5b=0.1d0
  k_hfc=0.0008d0
  k_init = 1.12d0 
  k_rids=10.0d0
  k_repl=1.12d0
  
  k_sp_5b=0.1d0
  k_sinit = 1.12d0 
  k_sids=10.0d0
  k_srepl=1.12d0

  !== assembly and other rates == !
  
  k_rip=0.6d0
  k_assembly=1.2d-7 !6.7d-10
  k_out_rp=0.307d0
  
  k_sip=0.6d0
  k_out_sp=0.307d0

  !==degradation rates == !

  k_deg_r=0.26d0
  k_deg_syn=0.26d0

!  k_deg_ns5b(:) = 0.0d0
  k_deg_si =0.10d0
  k_deg_ns=0.110d0
  k_deg_vms=0.001d0

  ! == populations == ! 

  rna(:)=0
  rna(1)=moi_v ! --  one infecting particle

  syn(:)=0

  ns5b(:)=0
  
  ribo = 5000
  tc = 0
  pp = 0
  s = 180*moi_v
  ns34a = 0
  ns5a = 0
  hf = 30
  hfc = 0 
 
  !== init ==! 

  b(:) = 0.0d0

  btot = 0.0d0

  time = 0.0d0

  time_end = 2400.5d0
  
  tout = 1d-6

  !iseed = 89693
  
  !read(*,*) iseed 
  !call srand(iseed) 

  r=0.0d0

  DO WHILE (time < time_end)  ! == start run == !
  
  IF (time .LE. 21.0d0) THEN
      k_deg_s=0.61d0
  ELSE
      k_deg_s=0.1d0
  ENDIF
  
  IF (time .GE. syntime .AND. synentry == 0) THEN
      syn(1) = moi_s
      synentry = 1
      s = s + 180*moi_s
  ENDIF
  
  ! == calc probabilities == !

  b(:) = 0.0d0

  !translation
 
  b(1) =  k_tc*rna(1)*ribo  ! rna + ribo
  b(2) =  k_trans*tc ! translation
  b(3) =  k_cleavage*pp ! cleavage

  !replication 

  b(4) = k_rp_5b * ns5b(1) * rna(1) ! ns5b + rna  
  b(5) = k_hfc * ns5a * hf ! hf + ns5a 
  b(6) = k_rip * hfc * rna(6) ! hfc + ns5b + rna 
  b(7) = k_init * rna(3) ! 
  b(8) = k_rids * rna(7) * ns5b(2) ! 
  b(9) = k_repl * rna(4) ! 

  !assembly and other 
  if (s>=180) then
  b(10) = k_assembly * rna(5) * s !
  else
  b(10) = 0.0d0
  endif
  b(11) = k_out_rp * rna(5) ! rna exiting VMS

  ! degradation  
  b(12) = k_deg_r * rna(1) !
  b(13) = k_deg_s * s ! 
  b(14) = k_deg_ns * ns34a ! 
  b(15) = k_deg_ns * ns5a ! 
  b(16) = k_deg_ns * ns5b(1) ! 
  b(17) = k_deg_vms * rna(3) ! 
  b(18) = k_deg_vms * rna(7) ! 
  b(19) = k_deg_vms * ns5b(2) ! 
  b(20) = k_deg_vms * rna(4) ! 
  b(21) = k_deg_vms * rna(5) ! 

  !synth replication 

  b(22) = k_sp_5b * ns5b(1) * syn(1) ! ns5b + rna  
  b(23) = k_sip * hfc * syn(6) ! hfc + ns5b + rna 
  b(24) = br * k_sinit * syn(3) ! 
  b(25) = k_sids * syn(7) * ns5b(2) ! 
  b(26) = br * k_srepl * syn(4) ! 

  !synth assembly and other 
  if (s>=180) then
  b(27) = k_assembly * syn(5) * s * ba
  else
  b(27) = 0.0d0
  endif
  b(28) = k_out_sp * syn(5) ! rna exiting VMS

  !synth degradation  
  b(29) = k_deg_syn * syn(1) !
  b(30) = k_deg_vms * syn(3) ! 
  b(31) = k_deg_vms * syn(7) ! 
  b(32) = k_deg_vms * syn(4) ! 
  b(33) = k_deg_vms * syn(5) ! 
  
  btot = 0.0d0

  do i=1,33
    btot = btot+b(i)
  end do 
  
  IF((btot).eq.0.0)EXIT

  r=RANDOM(iseed)
  
  tau=DLOG(1.0d0/r)/btot

  time = time + tau 
  
  IF ( time >= tout ) THEN
      
      WRITE(9,*) time,rna(8)
          
      tout = tout + 24.0d0

  ENDIF
  
  r=RANDOM(iseed)

  r = r*btot 

  i=0
  btot =0.0d0
  
  DO WHILE (btot < r)
    i = i + 1
    btot = btot +b(i)
  ENDDO

  IF ( i == 1) THEN 
    ! rna + ribo 
    rna(1) = rna(1) - 1 
    rna(2) = rna(2) + 1
    ribo = ribo - 1 
    tc = tc + 1 

  ENDIF 
  IF ( i == 2) THEN 
    !translation
    tc = tc -1 
    pp = pp +1 
    ribo = ribo + 1
    rna(1) = rna(1) + 1
    rna(2) = rna(2) - 1
  ENDIF

  IF ( i == 3) THEN 
    ! cleavage 
    pp = pp - 1
    s = s + 1
    ns34a = ns34a + 1 
    ns5a = ns5a + 1 
    ns5b(1) = ns5b(1) + 1 

  ENDIF

  IF ( i == 4) THEN 
    ! do this 
    ns5b(1) = ns5b(1) - 1 
    rna(1) = rna(1) - 1
    rna(6) = rna(6) + 1
  ENDIF
  
  IF ( i == 5) THEN 
    ! do this 
    ns5a = ns5a - 1 
    hf = hf - 1
    hfc = hfc + 1 
  ENDIF
  IF ( i == 6) THEN
    ! do this 
    hfc = hfc - 1 
    rna(6) = rna(6) - 1 
    rna(3) = rna(3) + 1
  ENDIF
  IF ( i == 7) THEN 
    ! do this 
    rna(3) = rna(3) - 1
    rna(7) = rna(7) + 1
    ns5b(2) = ns5b(2) + 1
    hf = hf + 1
  ENDIF
  IF ( i == 8) THEN 
    ! do this 
    rna(7) = rna(7) - 1 
    ns5b(2) = ns5b(2) - 1
    rna(4) = rna(4) + 1 
  ENDIF
  IF ( i == 9) THEN 
    ! do this 
    rna(4) = rna(4) - 1
    rna(7) = rna(7) + 1
    ns5b(2) = ns5b(2) + 1
    rna(5) = rna(5) + 1
  ENDIF
  IF ( i == 10) THEN 
    ! do this 
!    WRITE(*,*)'packaging' 
    rna(5) = rna(5) - 1 
    s = s - 180 
    rna(8) = rna(8) + 1
  !write(323,*)time,rna,ribo,tc,ns34a,ns5a,ns5b,s,syn,iseed
    ENDIF
  IF ( i == 11) THEN 
    ! do this 
    rna(5) = rna(5) - 1
    rna(1) = rna(1) + 1 
  ENDIF
  IF ( i == 12) THEN 
    ! do this 
    rna(1) = rna(1) - 1 
  ENDIF
  IF ( i == 13) THEN 
    ! do this 
    s = s -1 
  ENDIF
  IF ( i == 14) THEN 
    ! do this
    ns34a = ns34a - 1 
  ENDIF
  IF ( i == 15) THEN 
    ! do this 
    ns5a = ns5a - 1
  ENDIF
  IF ( i == 16) THEN 
    ! do this 
    ns5b(1) = ns5b(1) - 1 
  ENDIF
  IF ( i == 17) THEN 
    ! do this 
    rna(3) = rna(3) - 1
  ENDIF
  IF ( i == 18) THEN 
    ! do this 
    rna(7) = rna(7) - 1
  ENDIF
  IF ( i == 19) THEN 
    ! do this 
    ns5b(2) = ns5b(2) - 1
  ENDIF
  IF ( i == 20) THEN 
    ! do this 
    rna(4) = rna(4) - 1
  ENDIF
  IF ( i == 21) THEN 
    ! do this 
    rna(5) = rna(5) - 1 
  ENDIF
  IF ( i == 22) THEN 
    ! do this 
    ns5b(1) = ns5b(1) - 1 
    syn(1) = syn(1) - 1
    syn(6) = syn(6) + 1
  ENDIF
  IF ( i ==23) THEN
    ! do this 
    hfc = hfc - 1 
    syn(6) = syn(6) - 1 
    syn(3) = syn(3) + 1
  ENDIF
  IF ( i ==24) THEN 
    ! do this 
    syn(3) = syn(3) - 1
    syn(7) = syn(7) + 1
    ns5b(2) = ns5b(2) + 1
    hf = hf + 1
  ENDIF
  IF ( i ==25) THEN 
    ! do this 
    syn(7) = syn(7) - 1 
    ns5b(2) = ns5b(2) - 1
    syn(4) = syn(4) + 1 
  ENDIF
  IF ( i ==26) THEN 
    ! do this 
    syn(4) = syn(4) - 1
    syn(7) = syn(7) + 1
    ns5b(2) = ns5b(2) + 1
    syn(5) = syn(5) + 1
  ENDIF
  IF ( i == 27) THEN 
    ! do this 
    !WRITE(*,*)'packaging' 
    syn(5) = syn(5) - 1 
    s = s - 180 
    syn(8) = syn(8) + 1
  !write(323,*)time,rna,ribo,tc,ns34a,ns5a,ns5b,s,syn,iseed
    ENDIF
  IF ( i == 28) THEN 
    ! do this 
    syn(5) = syn(5) - 1
    syn(1) = syn(1) + 1 
  ENDIF
  IF ( i == 29) THEN 
    ! do this 
    syn(1) = syn(1) - 1 
  ENDIF
  IF ( i == 30) THEN 
    ! do this 
    syn(3) = syn(3) - 1
  ENDIF
  IF ( i == 31) THEN 
    ! do this 
    syn(7) = syn(7) - 1
  ENDIF
  IF ( i == 32) THEN 
    ! do this 
    syn(4) = syn(4) - 1
  ENDIF
  IF ( i == 33) THEN 
    ! do this 
    syn(5) = syn(5) - 1 
  ENDIF
  
  ENDDO 
  
  ENDDO

  !== write-out-final ==!

  !write(323,*)time,rna,ribo,tc,ns34a,ns5a,ns5b,s,syn,iseed
  END PROGRAM
  
  !!!!! Creating a Random Number !!!!!
      
      DOUBLE PRECISION FUNCTION RANDOM (ISEED)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(INOUT) :: iseed

        !=== VARIABLES ===!

        INTEGER :: hi,lo,test

        INTEGER, PARAMETER :: ae = 16807
        INTEGER, PARAMETER :: m = 2147483647
        INTEGER, PARAMETER :: q = 127773
        INTEGER, PARAMETER :: re = 2836


        hi = INT(iseed/q)
        lo = MODULO(iseed,q)

        test = ae * lo - re * hi

        IF ( test > 0 ) THEN

          iseed = test

        ELSE

          iseed = test + m

        ENDIF

        RANDOM = DBLE(iseed) / DBLE(m)

        RETURN

      END FUNCTION RANDOM
