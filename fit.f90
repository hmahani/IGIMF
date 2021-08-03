subroutine fit(alpha,slp,y_icp,slp_IGIMF,y_icp_IGIMF)
use utilmf
implicit none

integer::i,j
integer, intent(in)::alpha
real(8), intent(out)::slp,y_icp,slp_IGIMF,y_icp_IGIMF
real(8)::x_avg, y_avg, nslp,sigma_x,sigma_y,sigma_xy, sigma_x2, sigma_y_IGIMF, sigma_xy_IGIMF,y_avg_IGIMF
integer::j_d, j_u

	if (alpha .eq. 1)then
		j_d=5
		j_u=10
		nslp=6.0
	else if (alpha .eq. 2) then
		j_d=12
		j_u=14
		nslp=3.0
	else if (alpha .eq. 3) then
		j_d=13
		j_u=19
		nslp=7.0
	else if (alpha .eq. 4) then
		j_d=20
		j_u=25
		nslp=6.0
	end if


	sigma_x=0.0
	sigma_y=0.0
	sigma_xy=0.0
	sigma_x2=0.0
	sigma_y_IGIMF=0.0
	sigma_xy_IGIMF=0.0

	do i=j_d, j_u
		sigma_x=sigma_x+newlogm(i)
		sigma_y=sigma_y+log10(new_NIMF_f(i)/(newlogm(i)-newlogm(i-1)))
		sigma_y_IGIMF=sigma_y_IGIMF+log10(new_n_igimf_f(i)/(newlogm(i)-newlogm(i-1)))
	end do

	x_avg=sigma_x/nslp
	y_avg=sigma_y/nslp
	y_avg_IGIMF=sigma_y_IGIMF/nslp
		
	do i=j_d, j_u
		sigma_xy=sigma_xy+(newlogm(i)-x_avg)*(log10(new_NIMF_f(i)/(newlogm(i)-newlogm(i-1)))-y_avg)
		sigma_xy_IGIMF=sigma_xy_IGIMF+(newlogm(i)-x_avg)*log10(new_n_igimf_f(i)/(newlogm(i)-newlogm(i-1)))
		sigma_x2=sigma_x2+((newlogm(i)-x_avg)**2.0)
	end do

	slp=sigma_xy/sigma_x2
	y_icp=y_avg-slp*x_avg

	slp_IGIMF=sigma_xy_IGIMF/sigma_x2
	y_icp_IGIMF=y_avg_IGIMF-slp_IGIMF*x_avg

end subroutine
