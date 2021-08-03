subroutine binary(local_tho,iso)
use utilmf
implicit none

integer, intent(in)::iso,local_tho
 integer :: values(1:8), kk
 integer, dimension(:), allocatable :: seed
 real(8) :: r
 real(8)::l1, l2
 integer,parameter:: MyLong = selected_int_kind (12)
 integer (kind=MyLong) :: k, n_star,n_cnt, i,j,r1, r2, n_local


 call date_and_time(values=values)
 call random_seed(size=kk)
 allocate(seed(1:kk))
 seed(:) = values(8)
 call random_seed(put=seed)

	n_star=n_tot*1	
	n_cnt=(prcnt/100.0)*(n_tot/2.0)/1	!number of binaries

!write(*,*)iso,n_star,n_cnt
!write(*,*)

allocate (nb(0:n_t),nprim(0:n_t))
allocate ( total_lum_b(1:n_t))

	do j=0, n_t
		nb(j)=0.0
		nprim(j)=0.0
		ntmp(j)=n(j)
	end do
!__________________
!	BINARYMERGE
!__________________

	do j=1,n_t
		do i=1,25
			if (log10(m(j))> newlogm(i-1) .and. log10(m(j))<= newlogm(i)) then
				new_n(i)=new_n(i)+n(j)
				new_n_igimf(i)=new_n_igimf(i)+NIG(j)
				exit
			end if
		end do
	end do

	do i=1,25
		 new_NIMF_f(i)= new_NIMF_f(i)+new_n(i)
		 new_n_igimf_f(i)=new_n_igimf_f(i)+new_n_igimf(i)
	end do


DEALLOCATE (m)
DEALLOCATE (dm)
DEALLOCATE (dlm)
DEALLOCATE (n)
DEALLOCATE (phi)
DEALLOCATE (wght)
DEALLOCATE (wght_tmp)
DEALLOCATE (nb)
DEALLOCATE (nprim)
DEALLOCATE (ntmp)
DEALLOCATE (total_mass)
DEALLOCATE (total_lum)
DEALLOCATE (lumino)
DEALLOCATE (total_lum_b)
DEALLOCATE (phi_prim)
DEALLOCATE (PHI_IG)
DEALLOCATE (NIG)
DEALLOCATE (total_mass_IG)
DEALLOCATE (total_lum_IG)
DEALLOCATE (lU)
DEALLOCATE (lB)
DEALLOCATE (lV)
DEALLOCATE (lI)
DEALLOCATE (lJ)
DEALLOCATE (lK)
DEALLOCATE (lR)



end subroutine
