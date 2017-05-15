pro doublePendulum, _EXTRA = extra
  compile_opt idl2
  ireset, /NO_PROMPT
  
  outfile = getenv('user') + path_sep() + 'pendulums.gif'
  
  
  file_delete, outfile, /QUIET
  
  ;initial conditions
  th01 = 90  ;degrees
  dth01 = 0  ;degrees/second
  th02 = 0   ;degrees
  dth02 = -360  ;degrees/second
  
  ;pendulum parameters
  m1 = 0.5
  L1 = 0.5
  m2 = 2
  L2 = 0.5
  
  ;simlation parameters
  maxT = 10.0
  dT = 0.0001d
  n = ceil(double(maxT)/dT)
  g = 3d;9.81d
  
  ;animation parameters
  fps = 60
  dTGoal = 1d/30
  frameSpace = floor(dTGoal/dT) > 1
  waitTime = 1d/fps
  pathPrint = 5
  addPath = 0
  
  ;preallocate some arrays
  th1   = make_array(n, TYPE = 5, /NOZERO) 
  thd1  = make_array(n, TYPE = 5, /NOZERO)
  thdd1 = make_array(n, TYPE = 5, /NOZERO)
  th2   = make_array(n, TYPE = 5, /NOZERO)
  thd2  = make_array(n, TYPE = 5, /NOZERO)
  thdd2 = make_array(n, TYPE = 5, /NOZERO)
  Ttot  = make_array(n, TYPE = 5, /NOZERO)
  T1    = make_array(n, TYPE = 5, /NOZERO)
  T2    = make_array(n, TYPE = 5, /NOZERO)
  Utot  = make_array(n, TYPE = 5, /NOZERO)
  time  = findgen(n)*dT
  
  ;initialize our data
  th1[0]  = th01*!DTOR
  thd1[0] = dth01*!DTOR
  th2[0]  = th02*!DTOR
  thd2[0] = dth02*!DTOR
  
  ;loop over it all
  for i=0,n-1 do begin
    ;calculate our energies
    Ttot[i] = 0.5*m1*((thd1[i])^2)*(L1^2) + $
      0.5*m2*((thd1[i]^2)*(L1^2) + (thd2[i]^2)*(L2^2) + 2*thd1[i]*thd2[i]*L1*L2*cos(th1[i] - th2[i]))
    Utot[i] = -g*L1*(m1 + m2)*cos(th1[i]) - m2*L2*g*cos(th2[i])
    
    ; get individual rotational kinetic energies about their
    ; reference frames
;    T1[i] = .5*m1*(L1*thd1[i])^2 + .5*m1
;    T2[i] = .5*m2*(L2*thd2[i])^2
    
    ;calculate our first angular acceleration
    top1 = -m2*L1*(thd1[i]^2)*sin(th1[i] - th2[i])*cos(th1[i] - th2[i]) + $
      g*m2*sin(th2[i])*cos(th1[i] - th2[i]) - $
      m2*L2*(thd2[i]^2)*sin(th1[i] - th2[i]) - $
      (m1 + m2)*g*sin(th1[i])
    bottom1 = L1*(m1 + m2) - m2*L1*(cos(th1[i] - th2[i])^2)
    thdd1[i] = top1/bottom1
    
    ;calcualte our second angular acceleration
    top2 = m2*L2*(thd2[i]^2)*sin(th1[i] - th2[i])*cos(th1[i] - th2[i]) + $
      g*sin(th1[i])*cos(th1[i] - th2[i])*(m1 + m2) + $
      L1*(thd1[i]^2)*sin(th1[i] - th2[i])*(m1 + m2) - $
      g*sin(th2[i])*(m1 + m2)
    bottom2 = L2*(m1 + m2) - m2*L2*(cos(th1[i] - th2[i])^2)
    thdd2[i] = top2/bottom2
    
    ;make sure we stop at the right time so we don't illegally
    ;subscript our arrays
    if (i eq (n - 1)) then continue
    
    ;do our numeric integration
    thd1[i+1] = dT*thdd1[i] + thd1[i]
    thd2[i+1] = dT*thdd2[i] + thd2[i]
    th1[i+1]  = th1[i] + dT*thd1[i] + .5*(dT^2)*thdd1[i]
    th2[i+1]  = th2[i] + dT*thd2[i]  + .5*(dT^2)*thdd2[i]
  endfor
  
  ;convert to X and Y
  x1 = L1*sin(th1)
  y1 = -L1*cos(th1)
  x2 = x1 + L2*sin(th2)
  y2 = y1 - L2*cos(th2)
  
  ;get the total energy
  Etot = Ttot + Utot
  Eerror = Etot - Etot[0]

  ;make our plot
  range = (L1 + L2)*[-1,1]
  p = plot(x1, y1, COLOR = 'blue', XRANGE = range, YRANGE = range, ASPECT_RATIO = 1, /NODATA)
  dims = p.window.dimensions
  p1 = plot(x2, y2, XRANGE = range, YRANGE = range, CURRENT = p, /OVERPLOT, /NODATA)
  
  ;add positions of the balls
  line = plot([0, x1[0], x2[0]], [0, y1[0], y2[0]], CURRENT = p, /OVERPLOT, THICK = 2)
  b1 = bubbleplot(x1[0], y1[0], COLOR = 'blue', CURRENT = p, /OVERPLOT, SIZING = .2)
  b2 = bubbleplot(x2[0], y2[0], COLOR = 'black', CURRENT = p, /OVERPLOT, SIZING = .2)
  
  ;set up a color table!
  levels = 256.
  ctable = COLORTABLE([[0,0,255],[255,0,0]], NCOLORS = levels, /TRANSPOSE)
  
  ;animate our plot
  frames = [0:n-1:frameSpace]
  foreach i, frames, idx do begin
    ;get start of draw time
    tstart = systime(/SECONDS)
    
    ;stop updating the window
    p.window.refresh, /disable
    
    ;add path
    if (addPath) then begin
      if ((i gt 0) AND ~(idx mod pathPrint)) then begin
        pi1 = plot(x1[frames[idx-pathPrint]:i], y1[frames[idx-pathPrint]:i], COLOR = 'blue', CURRENT = p, /OVERPLOT)
        pi2 = plot(x2[frames[idx-pathPrint]:i], y2[frames[idx-pathPrint]:i], CURRENT = p, /OVERPLOT)
      endif
    endif
    
    ;add lines and pendulumns
    line.SetData, [0, x1[i], x2[i]], [0, y1[i], y2[i]]
    b1.SetData, x1[i], y1[i]
    b2.SetData, x2[i], y2[i]
    
;    if (color_set) then begin
;      b1.COLOR = ctable[*,body_colors]
;      b2.COLOR = ctable[*,body_colors]
;    endif
    
    ;update window
    p.window.refresh
    
    ;save our image frame
    if ~(idx mod 2) then begin
      img = p.window.copyWindow()
      ;convert the image to 8 bit
      img = color_quan(img, 1, r, g, b)
      ;img = congrid(img, dims[0]/2, dims[1]/2)
      ;duplicate frames otherwise low frame-rate can be a problem with some players
      WRITE_GIF, outfile, img, r, g, b, delay_time = 2.0*dTGoal*100 ,$
        /MULTIPLE, REPEAT_COUNT = 0
    endif

    
    ;get and of processing time and then wait what we need to
    tTot = systime(/SECONDS) - tStart
    if (tTot lt waitTime) then begin
      wait, (waitTime - tTot)
    endif
  endforeach

  ; Close the file and window
  WRITE_GIF, outfile, /CLOSE

end


doublePendulum, THIS = 1, THAT = 0

end
