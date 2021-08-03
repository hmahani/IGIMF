subroutine igimf(iso,i)
use utilmf
implicit none

integer::j,star_cnt, f
integer, intent(in)::i,iso
! integer,parameter:: MyLong = selected_int_kind (12)
! integer (kind=MyLong) :: k, n_star

      M_1=0.08
      M_2=0.5
      M_3=1.0
      M_4=150.0   
      a1=-1.3
      a2=-2.3
      a_1=-1.30-0.5*FeH
      a_2=-2.30-0.5*FeH


		

!													1

          Do f=1,n_igimf
                Mecl(f)=0.0
                Mecl(f)=10.0**(delta*f+log10(Mecl(0)))
        
              if(Mecl(f)<=Mecl_max) then

!													2	

                  Meclave(f)=(Mecl(f)+Mecl(f-1))/2.0
                  xx(f)=0.99*(0.61*log10(Mecl(f))+2.85-6.0)-0.14*FeH
                  dMecl(f)=Mecl(f)-Mecl(f-1) 
                  m_max(f)=10.0**(2.56*log10(Mecl(f))*(3.82**9.17+(log10(Mecl(f)))**9.17)**(-1.0/9.17)-0.38)
                  !write(*,*)Mecl(f),dMecl(f)
                  if(m_max(f)>150)then
                      m_max(f)=150
                  endif  
!													3

                  if (xx(f)>=-0.87 )then
                      a_3=(-1.0*(1.94-0.41*xx(f)))
                  else if(xx(f)<-0.87)then
                      a_3=-2.3
                  endif

!													4

                     k2= 1./(( 2.*(M_2**(a_1+2.)-M_1**(a_1+2.))/(a_1+2.)+ (M_3**(a_2+2.)-M_2**(a_2+2.))/(a_2+2.) &
                         &+ (m_max(f)**(a_3+2.)-M_3**(a_3+2.))/(a_3+2.)  ))
                     k1= 2.0*k2  
                     k3=k2  !since M3=1

!													5

                  if (m(i)>=0.08 .and. m(i)<0.5) then
                    phi_prim(i)=(Mecl(f)*K1)*m(i)**(a_1)
                  else if (m(i)>=0.5 .and. m(i)<1.0) then 
                    phi_prim(i)=(Mecl(f)*K2)*m(i)**(a_2)  
                  else if (m(i)>=1. .and. m(i)<150.0) then  
                    phi_prim(i)=(Mecl(f)*K3)*m(i)**(a_3)  
                  endif 
             
!													6

                  if(Bt==-2.)then
                   C1=1./(log(Mecl_max)-log(Mecl(0)))
                  else if(Bt>-2. .or. Bt<-2.)then 

                   C1=1./ ( ( Mecl_max**(Bt+2.)-Mecl(0)**(Bt+2.) )/(Bt+2.))
                  endif              
!													7
               
                  if(m(i)<=m_max(f)) then
                     PHI_IG(i)=phi_prim(i)*C1*Meclave(f)**(Bt)*dMecl(f)+PHI_IG(i)
                  elseif (m(i)>m_max(f)) then 
                     PHI_IG(i)=0.0
                  end if
        
             endif

          End do
!													8


             NIG(i)=PHI_IG(i)*dm(i)*SFR(iso)*(ddt(iso))

		lBIG_t(iso)=lBIG_t(iso)+lB(i)*NIG(i)
		lVIG_t(iso)=lVIG_t(iso)+lV(i)*NIG(i)
		lRIG_t(iso)=lRIG_t(iso)+lR(i)*NIG(i)
		lKIG_t(iso)=lKIG_t(iso)+lK(i)*NIG(i)
		lIIG_t(iso)=lIIG_t(iso)+lI(i)*NIG(i)
		lJIG_t(iso)=lJIG_t(iso)+lJ(i)*NIG(i)
		lUIG_t(iso)=lUIG_t(iso)+lU(i)*NIG(i)
                  
  
end subroutine
