      subroutine deltaa(n,dxi)

      include 'inclm1.f'

      real*8 dxi(*)
      real*8 dltau(mn), dlp(mn), dlu(mn), dlar(mn), dlg(mn), dlg1(mn)

      common /delta/ dltau, dlp, dlu, dlar, dlg, dlg1
      do 10 i = 2, n-1
         dxi1 = dxi(i+1) + dxi(i)
         dxim = dxi(i) + dxi(i-1)
         dxi2m = 2.d0*dxi(i-1) + dxi(i)
         dxi21 = 2.d0*dxi(i+1) + dxi(i)
         dxis = dxi(i-1)+dxi(i)+dxi(i+1)
         dpj = p(i+1) - p(i)
         dpjm = p(i) - p(i-1)
         if(dpj*dpjm.gt.0.d0)then
            dpmin2 = 2.d0*min(abs(dpj),abs(dpjm))
            dlp0 = dxi(i)/dxis*(dxi2m/dxi1*dpj+dxi21*dpjm/dxim)
            adlp0 = abs(dlp0)
            dlp(i) = min(adlp0,dpmin2)*dlp0/adlp0
         else
            dlp(i) = 0.d0
         endif
         duj = u(i+1) - u(i)
         dujm = u(i) - u(i-1)
         if(duj*dujm.gt.0.d0)then
            dumin2 = 2.d0*min(abs(duj),abs(dujm))
            dlu0 = dxi(i)/dxis*(dxi2m/dxi1*duj+dxi21*dujm/dxim)
            adlu0 = abs(dlu0)
            dlu(i) = min(adlu0,dumin2)*dlu0/adlu0
         else
            dlu(i) = 0.d0
         endif
         dtauj = tau(i+1) - tau(i)
         dtaujm = tau(i) - tau(i-1)
         if(dtauj*dtaujm.gt.0.d0)then
            drmin2 = 2.d0*min(abs(dtauj),abs(dtaujm))
            dltu0 = dxi(i)/dxis*(dxi2m/dxi1*dtauj+dxi21*dtaujm/dxim)
            adltu0 = abs(dltu0)
            dltau(i) = min(adltu0,drmin2)*dltu0/adltu0
         else
            dltau(i) = 0.d0
         endif
         darj = ar(i+1) - ar(i)
         darjm = ar(i) - ar(i-1)
         if(darj*darjm.gt.0.d0)then
            drmin2 = 2.d0*min(abs(darj),abs(darjm))
            dlar0 = dxi(i)/dxis*(dxi2m/dxi1*darj+dxi21*darjm/dxim)
            adlar0 = abs(dlar0)
            dlar(i) = min(adlar0,drmin2)*dlar0/adlar0
         else
            dlar(i) = 0.d0
         endif
         dgj = g(i+1) - g(i)
         dgjm = g(i) - g(i-1)
         if(dgj*dgjm.gt.0.d0)then
            dgmin2 = 2.d0*min(abs(dgj),abs(dgjm))
            dlg0 = dxi(i)/dxis*(dxi2m/dxi1*dgj+dxi21*dgjm/dxim)
            adlg0 = abs(dlg0)
            dlg(i) = min(adlg0,dgmin2)*dlg0/adlg0
         else
            dlg(i) = 0.d0
         endif
         dg1j = g1(i+1) - g1(i)
         dg1jm = g1(i) - g1(i-1)
         if(dg1j*dg1jm.gt.0.d0)then
            dg1mn2 = 2.d0*min(abs(dg1j),abs(dg1jm))
            dlg10 = dxi(i)/dxis*(dxi2m/dxi1*dg1j+dxi21*dg1jm/dxim)
            adlg10 = abs(dlg10)
            dlg1(i) = min(adlg10,dg1mn2)*dlg10/adlg10
         else
            dlg1(i) = 0.d0
         endif
 10   continue
      do 20 i = 1, 2
         dlp(i) = -dlp(5-i)
         dlar(i) = dlar(5-i)
         dlu(i) = dlu(5-i)
         dltau(i) = -dltau(5-i)
         dlg(i) = -dlg(5-i)
         dlg1(i) = -dlg1(5-i)
 20   continue
      dlp(n) = dlp(n-2)
      dlar(n) = dlar(n-2)
      dlu(n) = dlu(n-2)
      dltau(n) = dltau(n-2)
      dlg(n) = dlg(n-2)
      dlg1(n) = dlg1(n-2)
      return
      end
