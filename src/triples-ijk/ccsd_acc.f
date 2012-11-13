      subroutine terms_32_33_34_35(o1, o2, t2, vm, rl, tmp)
      use common_cc
      implicit none
      double precision o1(nou), o2(no2u2), t2(no2u2), vm(no3u),
     &                 rl(no2u2), tmp(nutr*nu)
      double precision, allocatable :: ve(:)
      integer :: i,a,j,jlo,jhi,nr,sr,ioff,iof1,iof2

      if(nutr*nu .lt. nou2) then
         write(6,*) 'tmp is sized incorrectly'
      endif
   
      allocate( ve(1:nu3) )

!$acc data present( o1, o2, t2, vm, rl )
!$acc&     create( tmp(1:nutr*nu), ve(1:nu3) )

!$acc parallel loop private(i) async(1)
      do i = 1,no2u2
         rl(i) = zero
      end do
!$acc end parallel loop

      CALL TRANSQ_ACC(T2,NOU)
      CALL INSITU_ACC(NO,NU,NU,NO,tmp,O2,12)
      CALL INSITU_ACC(NU,NO,NU,NO,tmp,O2,23)
      CALL TRANMD_ACC(VM,NO,NO,NO,NU,312)

      if(ddi_me.ne.0) then
!$acc parallel loop private(i) async(1)
         do i = 1,no3u
            vm(i) = zero
         end do
!$acc end parallel loop
      endif

      CALL DIV_EVEN(NO,DDI_NN,DDI_MY,NR,SR)
      JLO = SR - 1
      JHI = SR + NR - 2
 
      do j = jlo, jhi

         ioff = nou2*j + 1

! get and xfer ve to the gpu
         call ddi_get(d_vvvo,1,nutr,j*nu,((J*NU+NU-1)),tmp)
!$acc update device( tmp(1:nutr*nu) ) async(1)
         call ddcc_t_getve_acc(1,nu,j,tmp,ve)
!$acc wait(1)

!$acc host_data use_device( o1, ve, rl )
         CALL DGEMM(
     &             'N','N',NO,NU2,NU,ONE,O1,NO,ve,NU,
     &             ZERO,RL(IOFF),NO)
!$acc end host_data 
!$acc update host( rl(ioff:ioff+nou2-1) ) async(2)

         CALL TRANMD_ACC(ve,NU,NU,NU,1,13)

         iof1 = j*no2 + 1
!$acc host_data use_device( o2, ve, tmp )
         CALL DGEMM(
     &             'T','N',NO2,NU,NU2,ONE,O2,NU2,ve,NU2,
     &             ZERO,tmp,NO2)
!$acc end host_data 

         iof1 = j*no2 + 1
         iof2 = 1
!$acc parallel loop private(a,i)
         do a = 1,nu
            do i = 1,no2
               vm(iof1+i) = vm(iof1+i) + tmp(iof2+i)
            end do
            iof1 = iof1 + no3
            iof2 = iof2 + no2
         end do
!$acc end parallel loop
      end do

!$acc update host( vm(1:no3u) ) async(1)

!$acc wait(2)
      CALL DDI_MASTERS_GSUMF(305,RL,NO2U2)
!$acc update device( rl(1:no2u2) ) async(2)

!$acc wait(1)
      CALL DDI_MASTERS_GSUMF(306,VM,NO3U)
!$acc update device( vm(1:no3u) ) async(1)

!$acc wait(2)
      call transq_acc(rl,nou)
!$acc update host( rl(1:no2u2) ) async(99)
!$acc wait(1)
      call tranmd_acc(vm,no,no,no,nu,23)
!$acc end data
      deallocate( ve )
      end subroutine terms_32_33_34_35
