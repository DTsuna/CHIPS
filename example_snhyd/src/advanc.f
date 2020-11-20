      subroutine advanc(n,alpha,nadd,dt,dmass,encmg,time,
     $               boundr,e_charge_tot,injection_time,innerCell)

      include 'inclm1.f'
      include 'inclold.f'
      parameter ( pi = 3.141592654 )
      real*8 dmass(*), encmg(*)
      real*8 ps(mn), us(mn)
      real*8 nu(mn)
      real*8 vis, odt

      real*8 time,boundr
      real*8 amp
      real*8 period
      integer num_of_vib

      real*8 e_charge_tot
      real*8 e_charge_now
      real*8 e_charge_this
      integer e_in_cell
      real*8 injection_time
      integer kk
      integer innerCell
      common /riem/ ps, us
      data vis, odt / 0.1d0, 0.d0 /
      ial1 = int(alpha+1.d0)
      call old(n, nu)
      nu(2) = vis*max(u(2)-u(3),0.d0)


      e_in_cell = 10

      write(*,*)"innerCell in advanc.f =",innerCell

!     rad(2) = orad(2)+dt*(us(2)+nu(2))-odt*onu(2)
      rad(2) = 1d-15

      rad(innerCell) = boundr

c$$$      rad(2) = orad(2)+dt*us(2)
      do j = innerCell+1, n 		!from 3 to 4
         nu(j) = vis*max(u(j)-u(j+1),0.d0)
         rad(j) = orad(j)+dt*(us(j)+nu(j))-odt*onu(j)
         radm = 0.5*(rad(j)+rad(j-1))
         drad = rad(j)-rad(j-1)
         test = drad/radm
         if(test.lt.1d-5)then
            aa1 = (alpha+1)*radm**alpha*drad
         else
            a = rad(j)**int(alpha+1.d0)
            a1 = rad(j-1)**int(alpha+1.d0)
            aa1 = a-a1
         end if
         tau(j) = 4.d0*pi*aa1/((alpha+1.d0)*dmass(j))
      enddo
      tau(1) = tau(4)
      tau(2) = tau(3)
c$$$      tau(n-1) = tau(n-2)
c$$$      tau(n) = tau(n-2)
      rad(1) = -rad(3)
c$$$      rad(n) = orad(n)+dt*us(n)
c$$$      a = rad(n)**ial1
c$$$      a1 = rad(n-1)**ial1
c$$$      aa1 = a-a1
c$$$      tau(n) = 4.d0*pi*aa1/((alpha+1.d0)*dmass(n))
      
      call grav(n,encmg)
      do 20 j = innerCell+1, n	!from 3 to 4
         r2 = rad(j)*rad(j)
         r21 = rad(j-1)*rad(j-1)
         ar(j) = 4.*pi*(r2+r21+rad(j)*rad(j-1))/3.d0
 20   continue

      u(innerCell)  = 0.d0
      us(innerCell) = 0.d0
!      ps(3) = (p(3)+p(4))/2	!riemntでのps(3)はそのまま使えないので
      !ps(3) = 1.097d+23	!riemntでのps(3)はそのまま使えないので
!mesaでの値を時間変化を含めて

      do 30 j = innerCell+1, n	!from 3 to 4
         r2 = rad(j)*rad(j)
         or2 = orad(j)*orad(j)
         r21 = rad(j-1)*rad(j-1)
         or21 = orad(j-1)*orad(j-1)
         ab = (r2+or2+rad(j)*orad(j))/3.d0
         ab1 = (r21+or21+rad(j-1)*orad(j-1))/3.d0
         grvm = 0.5d0*(grv(j)+ogrv(j))
c$$$         u(j) = ou(j)+2.*pi*(ab+ab1)*dt/dmass(j)
c$$$     $        *(ps(j-1)+nu(j-1)*ou(j-1)/otau(j-1)
c$$$     $        -ps(j)-nu(j)*ou(j)/otau(j))
c$$$     $        + dt*grvm 
c$$$         e(j) = oe(j)+4.*pi*dt/dmass(j)*(ab1*(us(j-1)*ps(j-1)
c$$$     $        +nu(j-1)*oe(j-1)/otau(j-1))
c$$$     &        -ab*(us(j)*ps(j)+nu(j)*oe(j)/otau(j)))
         u(j) = ou(j)+2.*pi*(ab+ab1)*dt/dmass(j)
     $        *(ps(j-1)-ps(j))
     $        +dt*grvm
         e(j) = oe(j)+4.*pi*dt/dmass(j)*(ab1*us(j-1)*ps(j-1)
     $        -ab*us(j)*ps(j))

         if(j.eq.10)then
           write(*,*)"====== advanc report at 10 ====="
           write(*,*)"u in first=", 2.*pi*(ab+ab1)*dt/dmass(j)*
     $           (ps(j-1)-ps(j))
           write(*,*)"u in second=",dt*grvm
         end if
 30   continue

      if(time.lt.5.d1.and.time+dt.gt.5.d1)then
        e_charge = (e_charge_tot/injection_time)*(time+dt-5.d1)
        write(*,*)"e_charge=",e_charge
        do kk = 10,9+e_in_cell
          e(kk) = e(kk) + (e_charge/e_in_cell)/dmass(kk)
        end do
      end if 


      if(time.gt.5.d1.and.time.lt.5.d1+injection_time)then
        e_charge = (e_charge_tot/injection_time)*dt
        write(*,*)"e_charge=",e_charge
        if(time+dt.gt.5.d1+injection_time)then
          e_charge = (e_charge_tot/injection_time)*
     $    (5.d1+injection_time-time)
          write(*,*)"e_charge=",e_charge
        end if
        do kk = 10,9+e_in_cell 
          e(kk) = e(kk) + (e_charge/e_in_cell)/dmass(kk)
        end do
      end if




      !write(*,*)'u(4)e(4)'
      !print *,ps(3)
      !print *,ps(4)

c$$$      u(n) = u(n-1)
c$$$      e(n) = e(n-1)

      do 40 i = 1, 2
         u(i) = -u(5-i)
         ar(i) = -ar(5-i)
         e(i) = e(5-i)
c         u(i+n-2) = u(n-2)
c         e(i+n-2) = e(n-2)
 40   continue
      odt = dt
c$$$      do j = n-nadd, n-1
c$$$         test = 2.d0*tau(j)/(tau(j-1)+tau(j+1))
c$$$         if(test.lt.1d-1)then
c$$$            rad(j) = (0.5*rad(j)+0.5*rad(j+1))
c$$$            rad(j-1) = (0.5*rad(j-1)+0.5*rad(j-2))
c$$$         end if
c$$$      enddo
      return
      end
