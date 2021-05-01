      subroutine view(nn,idev,time,rad,tau,p,u,ye,lum,tmp)
      integer n, idev, init, fileplot
      real*8 rad(*), tau(*), p(*), u(*), ye(*), lum(*), tmp(*), solm,
     $     time
      include 'inclmn.f'
      parameter ( solm = 1.989d33 )
      real gm(mn), grho(mn),gp(mn),gv(mn),gt(mn),gy(mn),gl(mn),
     $     grho_min, grho_max, gp_min, gp_max, gv_min, gv_max, 
     $     gt_min, gt_max, gl_min, gl_max, gy_min, gy_max ,
     $     xl1, xl2, yl1 , yl2
      character*80 ctime
      data init /1/
      write(ctime,'(''\\r t = '',1pe11.4,'' s\\e'')')time
      if(init.eq.1)then
         call device(idev)
         call tsetup
         if(idev.eq.11)call setphysical(0.,0.,639.,741.)
         call getloc(xl1,yl1,xl2,yl2)
c         write(*,'(1p4e12.4)')xl1,yl1,xl2,yl2*1.5
         call setloc(xl1,yl1,xl2,yl2*1.5)
         init = 0
         return
      end if
      do 10 j = 1, nn-2
         gm(j) = log10(abs(rad(j+2)))
         grho(j) = -log10(tau(j+2))
         gp(j) = log10(p(j+2))
         gv(j) = u(j+2)*1d-8
         gy(j) = ye(j+2)
c         gl(j) = lum(j+2)*1d-40
         gl(j) = u(j+2)
         gt(j) = log10(tmp(j+2))
 10   continue
      n = nn - 2
      grho_min = 1e30
      grho_max = -1e30
      gp_min = 1e30
      gp_max = 0
      gv_min = 1e30
      gv_max = 0
      gy_min = 1e30
      gy_max = 0
      gt_min = 1e30
      gt_max = 0
      gl_min = 1e30
      gl_max = -1e30
      do 20 j = 1, n-2
         grho_min = min(grho_min,grho(j))
         grho_max = max(grho_max,grho(j))
         gp_min = min(gp_min,gp(j))
         gp_max = max(gp_max,gp(j))
         gv_min = min(gv_min,gv(j))
         gv_max = max(gv_max,gv(j))
         gt_min = min(gt_min,gt(j))
         gt_max = max(gt_max,gt(j))
         gy_min = min(gy_min,gy(j))
         gy_max = max(gy_max,gy(j))
         gl_min = min(gl_min,gl(j))
         gl_max = max(gl_max,gl(j))
 20   continue
      grho_max = grho_max*0.9
      gp_max = gp_max*1.1
      gv_max = gv_max*1.1
      gt_max = gt_max*1.1
      gy_max = gy_max*1.1
      gl_max = gl_max*1.1
      grho_min = grho_min*1.1
      gp_min = gp_min*0.9
      gv_min = gv_min*0.9
      gt_min = gt_min*0.9
      gy_min = gy_min*0.9
      gl_min = gl_min*0.9
c      write(*,*)gl_min,gl_max
      call erase
C---- Set the user coordinates and make a coordinate box
      call submargins(1.0,0.0)
      call window(1,5,1)
      call setlim(gm(1),grho_min, gm(n-2)*1.01, grho_max)
      call box(1,2)
      call xlabel(80,'Log R\e')
      call ylabel(80,'log \gr \e')
      call points(173.03,1,gm,grho,mn)
c      call connect(gm, grho, n)
      call window(1,5,2)
      call setlim(gm(1),gp_min, gm(n-2)*1.01, gp_max)
      call box(0,2)
      call ylabel(80,'\\rlog P \e')
      call points(173.03,1,gm,gp,mn)
c      call connect(gm, gp, n)
      call window(1,5,3)
      call setlim(gm(1),gv_min, gm(n-2)*1.01, gv_max)
      call box(0,2)
      call ylabel(80,'\\rV\d8 \e')
      call points(173.03,1,gm,gv,mn)
c      call connect(gm, gv, n)
      call window(1,5,4)
      call setlim(gm(1),gt_min, gm(n-2)*1.01, gt_max)
      call box(0,2)
      call ylabel(80,'\\rLog T\e')
      call points(173.03,1,gm,gt,mn)
c      call connect(gm, gt, n)
      call window(1,5,5)
      call setlim(gm(1),gl_min, gm(n-2)*1.01, gl_max)
      call box(0,2)
      call ylabel(80,'\\rLum\\d40\e')
c      call connect(gm, gl, n)
      call points(173.03,1,gm,gl,mn)
      call tlabel(80,ctime)
      if(idev.eq.0) then
         write(*,'(i6)')fileplot(0)
      else
         call tidle
      end if
      return
      end
