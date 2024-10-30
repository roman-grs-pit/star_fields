import numpy as np
from starcount_NEW import get_starcounts as starcount
from astropy.table import Table

def insertstars(inputRA,inputDEC,seed=42,out_fn='teststars.txt'):
    area=0.281   #; Area in square degrees simulated. 0.016 is one field.
    pixscale=0.11 #   ; Pixel scale in arcseconds
    wavelength=1.58 #; Desired wavelength of output extinction, A_lambda, and
                    #;    must be input in microns.
    maglow=25.0     #; Faintest magnitude to be considered for stars
    maghi=5.0       #; Brightest magnitude to be considered for stars
    scale=1.0       #; 1 gives a square, 2.0 gives 1x2 rectangle, 3.0 gives 1x3 rectangle, etc.
    count=1         #; Starting number for output table
    template_starting_point=1   #; Row numbver of input SED list. For 28 stars at G0V (need to update)
    template_index=29 #; Number of absolute bins (each will be attached to a MS and Giant/WD if appropriate). Right now hard written in starcount pro
    extinct_bin=3 #; Number of bins extinction (Nominally at Hband) split into. A_H=[0,0.5,1.0]
    
    totalarea=area*3600.0*3600.0 #; Area in square arcseconds
    xdim=np.sqrt(scale)*np.sqrt(totalarea)/pixscale
    ydim=xdim/scale
    xdim2=xdim/2.0
    ydim2=ydim/2.0

    stardensity1,starmags = starcount(1,inputRA,inputDEC,wavelength)
    stardensity=area*stardensity1
    mag_delta=starmags[1]-starmags[2]
    sel = ((starmags >= maghi) & (starmags <= maglow))
    starmags2=starmags[sel]
    star_lth=len(starmags2)
    starcounts=np.zeros((template_index,extinct_bin,star_lth,2),dtype=int)
    stardensity = stardensity[:,:,sel,:] #select down to the magnitude range
    star_lth3 = 0 
    star_lth2 = 0 
    n0 = 0
    for j in range(0,template_index):
        for k in range(0,extinct_bin):
            for i in range(0,star_lth):         
                starnum=round(stardensity[j][k][i][0])
                if starnum == 0:
                    n0 += 1
                starcounts[j][k][i][0]=starnum
                star_lth2 += starnum
                starnum=round(stardensity[j][k][i][1])
                if starnum == 0:
                    n0 += 1
                
                starcounts[j][k][i][1]=starnum
                star_lth3 += starnum
#  ;  forprint,starmags(gg),starcounts(10,0,*,0),starcounts(10,0,*,1)
    
    #star_lth2=total(starcounts(*,*,*,0))
    #not sure of the quick way for python, just doing naive sum above
    #star_lth3=total(starcounts(*,*,*,1))
    print(star_lth2,star_lth3,star_lth2+star_lth3,n0)
    star_total=star_lth2+star_lth3
    star_input=[star_lth2,star_lth3]
    star_list=np.zeros((16,star_total))
    
    count=0
    rng = np.random.default_rng(seed=seed)
    check=rng.random(int(2*(star_total)))
    for j in range(0,template_index): 
        for k in range(0,extinct_bin):
            for l in range(0,2):
                star_index=int(j+template_index*l+(2*template_index)*k)# + template_starting_point)  
                #print(star_index,j,l,k,template_starting_point)
                for i in range(0,star_lth):
               
                    runningtot=0              
                    amount=starcounts[j][k][i][l]
               

                    while runningtot < amount:
                        x_new=xdim*check[count]
                        y_new=ydim*check[star_total+count]
                        newDEC=((y_new-ydim2)*pixscale/3600.0)+inputDEC
                        radDec=(newDEC/180.0)*np.pi
                        newRA=inputRA-(((x_new-xdim2)*pixscale/3600.0)/np.cos(radDec))
                        #print(star_index,j,l,k,template_starting_point)
                        star_list[0][count]=count
                        star_list[1][count]=x_new
                        star_list[2][count]=y_new
                        star_list[3][count]=8.0
                        star_list[4][count]=8.0
                        star_list[5][count]=0.0
                        star_list[6][count]=star_index
                        star_list[7][count]=0.0
                        star_list[8][count]=starmags2[i]
                        star_list[9][count]=0.0
                        star_list[10][count]=1
                        star_list[11][count]=0.0
                        star_list[12][count]=0.0
                        star_list[13][count]=star_list[8][count]
                        star_list[14][count]=newRA 
                        star_list[15][count]=newDEC
                        count=count+1
                        runningtot=runningtot+1
    startab = Table()
    names = ['index','Xpos','Ypos','8','8again','0','star_template_index','0again','magnitude','0again2','1','0again3','0again4','mag_again','RA','DEC']
    for i in range(0,len(names)):
        startab[names[i]] = star_list[i]
    startab.write(out_fn.replace('.txt','.ecsv'),overwrite=True)
    star_array = np.array(startab)
    h = ''
    for name in names:
        h += name+' '
    np.savetxt(out_fn,star_array,header=h)
    return True
