      subroutine advanc(n,alpha,nadd,dt,dmass,encmg,time,relt,
     $               boundr,e_charge_tot,injection_time,innerCell)

      include 'inclm1.f'
      include 'inclold.f'
      parameter ( pi = 3.141592654 )
      real*8 dmass(*), encmg(*)
      real*8 ps(mn), us(mn)
      real*8 nu(mn)
      real*8 vis, odt

      real*8 relt
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

      rad(2) = 1d-15

      rad(innerCell) = boundr

      do j = innerCell+1, n
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
      rad(1) = -rad(3)

      call grav(n,encmg)
      do 20 j = innerCell+1, n
         r2 = rad(j)*rad(j)
         r21 = rad(j-1)*rad(j-1)
         ar(j) = 4.*pi*(r2+r21+rad(j)*rad(j-1))/3.d0
 20   continue

      u(innerCell)  = 0.d0
      us(innerCell) = 0.d0
      !ps(3) = 1.097d+23

      do 30 j = innerCell+1, n
         r2 = rad(j)*rad(j)
         or2 = orad(j)*orad(j)
         r21 = rad(j-1)*rad(j-1)
         or21 = orad(j-1)*orad(j-1)
         ab = (r2+or2+rad(j)*orad(j))/3.d0
         ab1 = (r21+or21+rad(j-1)*orad(j-1))/3.d0
         grvm = 0.5d0*(grv(j)+ogrv(j))
         u(j) = ou(j)+2.*pi*(ab+ab1)*dt/dmass(j)
     $        *(ps(j-1)-ps(j))
     $        +dt*grvm
         e(j) = oe(j)+4.*pi*dt/dmass(j)*(ab1*us(j-1)*ps(j-1)
     $        -ab*us(j)*ps(j))

 30   continue

      if(time.lt.relt.and.time+dt.gt.relt)then
        e_charge = (e_charge_tot/injection_time)*(time+dt-relt)
        write(*,*)"e_charge=",e_charge
        do kk = 4,3+e_in_cell
          e(kk) = e(kk) + (e_charge/e_in_cell)/dmass(kk)
        end do
      end if


      if(time.gt.relt.and.time.lt.relt+injection_time)then
        e_charge = (e_charge_tot/injection_time)*dt
        write(*,*)"e_charge=",e_charge
        if(time+dt.gt.relt+injection_time)then
          e_charge = (e_charge_tot/injection_time)*
     $    (relt+injection_time-time)
          write(*,*)"e_charge=",e_charge
        end if
        do kk = 4,3+e_in_cell
          e(kk) = e(kk) + (e_charge/e_in_cell)/dmass(kk)
        end do
      end if


      do 40 i = 1, 2
         u(i) = -u(5-i)
         ar(i) = -ar(5-i)
         e(i) = e(5-i)
 40   continue
      odt = dt
      return
      end
