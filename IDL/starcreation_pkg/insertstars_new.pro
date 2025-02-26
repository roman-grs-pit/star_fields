pro insertstars,inputRA,inputDEC,output
    area=0.281   ; Area in square degrees simulated. 0.016 is one field.
    pixscale=0.11   ; Pixel scale in arcseconds
    wavelength=1.58 ; Desired wavelength of output extinction, A_lambda, and
                    ;    must be input in microns.
    maglow=25.0     ; Faintest magnitude to be considered for stars
    maghi=5.0       ; Brightest magnitude to be considered for stars
    scale=1.0       ; 1 gives a square, 2.0 gives 1x2 rectangle, 3.0 gives 1x3 rectangle, etc.
    count=1         ; Starting number for output table
    template_starting_point=1   ; Row numbver of input SED list. For 28 stars at G0V (need to update)
    template_index=29 ; Number of absolute bins (each will be attached to a MS and Giant/WD if appropriate). Right now hard written in starcount pro
    extinct_bin=3 ; Number of bins extinction (Nominally at Hband) split into. A_H=[0,0.5,1.0]
    
    totalarea=area*3600.0*3600.0 ; Area in square arcseconds
    xdim=sqrt(scale)*sqrt(totalarea)/pixscale
    ydim=xdim/scale
    xdim2=xdim/2.0
    ydim2=ydim/2.0

    starcount,inputRA,inputDEC,wavelength,1,starmags,stardensity1
    stardensity=area*stardensity1
    mag_delta=starmags(1)-starmags(2)
    gg=where(starmags ge maghi and starmags le maglow)
    starmags2=starmags(gg)
    star_lth=n_elements(gg)
    starcounts=intarr(template_index,extinct_bin,star_lth,2) 
    for j=0,template_index-1 do begin
       for k=0,extinct_bin-1 do begin
       for i=0,star_lth-1 do begin
          starnum=round(stardensity(j,k,gg(i),0))
          starcounts(j,k,i,0)=starnum
          starnum=round(stardensity(j,k,gg(i),1))
          starcounts(j,k,i,1)=starnum
       endfor
       endfor
    endfor
  ;  forprint,starmags(gg),starcounts(10,0,*,0),starcounts(10,0,*,1)
    star_lth2=total(starcounts(*,*,*,0))
    star_lth3=total(starcounts(*,*,*,1))
    print,star_lth2,star_lth3,star_lth2+star_lth3
    star_total=star_lth2+star_lth3
    star_input=[star_lth2,star_lth3]
    star_list=fltarr(16,star_total)
    
    count=long(0)
    seed = systime(1)
    check=randomu(seed,2*(star_total))
    for j=0,template_index-1 do begin
       for k=0,extinct_bin-1 do begin
          for l=0,1 do begin
             star_index=j+template_index*l+(2*template_index)*k + template_starting_point  
             
             for i=0,star_lth-1 do begin
               
               runningtot=0              
               amount=starcounts(j,k,i,l)
               

               while runningtot lt amount do begin
                  x_new=xdim*check(count)
                  y_new=ydim*check(star_total+count)
                  newDEC=((y_new-ydim2)*pixscale/3600.0)+inputDec
                  radDec=(newDec/180.0)*!pi
                  newRA=inputRA-(((x_new-xdim2)*pixscale/3600.0)/cos(radDec))
                  
                  star_list(0,count)=count
                  star_list(1,count)=x_new
                  star_list(2,count)=y_new
                  star_list(3,count)=8.0
                  star_list(4,count)=8.0
                  star_list(5,count)=0.0
                  star_list(6,count)=star_index
                  star_list(7,count)=0.0
                  star_list(8,count)=starmags2(i)
                  star_list(9,count)=0.0
                  star_list(10,count)=1
                  star_list(11,count)=0.0
                  star_list(12,count)=0.0
                  star_list(13,count)=star_list(8,count)
                  star_list(14,count)=newRA 
                  star_list(15,count)=newDec 
                  count=count+1
                  runningtot=runningtot+1
               endwhile
               
           endfor
         endfor
       endfor
    endfor
    
     writecol,output,long(star_list(0,*)),double(star_list(1,*)),double(star_list(2,*)),star_list(3,*),star_list(4,*),star_list(5,*),star_list(6,*),$
             star_list(7,*),star_list(8,*),star_list(9,*),star_list(10,*),star_list(11,*),star_list(12,*),$
             star_list(13,*),star_list(14,*),star_list(15,*)
end
