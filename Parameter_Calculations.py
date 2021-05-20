"""CODE AUTHOR: JORDAN SHERIDAN"""
import numpy as np
import pandas as pd
import math
import sys
import time
import json
import pprint
import pyvo as vo
from tkinter import messagebox, Toplevel, LabelFrame, Entry, Label, Button, IntVar

def MainProgram(Star,checked_stars,Results,runningpage):
    global period, period_err, Drop_in_flux, Drop_in_flux_error, data_frame 
    global exoarchive, outData_tic, outData_gaia, mode, second_attempt, M,r_star,v_star,L,M_err,r_star_err,v_star_err,L_err, star_variables
    mode = runningpage
    second_attempt = False
    data_frame = pd.DataFrame(columns = ['Star', 'Star Radius', 'Planet Classification', 'Period (d)','Period Error','Drop in flux','Drop in flux Error','Orbital radius (AU)','Orbital radius Error','Transit time (hrs)','Transit Error',
                                     'Planetary mass (M⊕)','Mass Error1','Mass Error2','Planetary density (g/cm\u00b3)','Density Error1','Density Error2','Planetary radius (R⊕)','Planetary radius Error',
                                     'Planetary velocity (km/s)','Velocity Error','Surface temperature (K)','Surface temperature Error','Surface gravity (m/s\u00b2)', 'Surface gravity Error1','Surface gravity Error2','Within habitable zone'])

    for df in Results:
        for i in range(len(df)):
            period = df['Period'][0]
            period_err = df['Period Error'][0]
            Drop_in_flux = df['Drop in Flux'][0]
            Drop_in_flux_error = df['Drop in Flux Error'][0]
            
            if period == 0:
                M,M_err,r_star,r_star_err,v_star,v_star_err,L,L_err = [None,None,None,None,None,None,None,None]
                star_variables = [None,None,None,None]
                NoResults(Star)
                return data_frame
            
            already_checked = [i for i in range(len(checked_stars)) if checked_stars[i][0] == Star]
            if len(already_checked) == 0:
                QueryParameters(Star)
                checked_stars.append([Star,outData_tic,outData_gaia,exoarchive])
                CalculateParameters(Star)
            else:
                outData_tic = checked_stars[already_checked[0]][1]
                outData_gaia = checked_stars[already_checked[0]][2]
                exoarchive = checked_stars[already_checked[0]][3]
                CalculateParameters(Star)
    
    return data_frame

#-------------------------------------------------------------------------#
#Jordan's code
#-------------------------------------------------------------------------#
#Significant figures
#-----------------------------------------------------------------------#
#ROUNDING FUNCTION
def significant_figures(value,error):
    if value == 0 or error == 0:
        sf_value = 0
        sf_error = 0
    else:
        if error == None:
            return value, None
        elif type(error) == list:
            sigfigs = min([-int(math.floor(math.log10(abs(err)))) for err in error])
            sf_error = [round(error[0],sigfigs),round(error[1],sigfigs)]
            sf_error = [int(10**-sigfigs) if x == 0 else x for x in sf_error]
            if sigfigs <= 0: sf_error = [int(err) if int(err) != 0 else int(10**-sigfigs) for err in sf_error]
        else:
            sigfigs = - int(math.floor(math.log10(abs(error))))
            sf_error = round(error,sigfigs)
            if sigfigs <= 0: sf_error = int(sf_error)
            
        sf_value = round(value,sigfigs)
        if sigfigs <= 0: sf_value = int(sf_value)
        if sf_value == 0: sf_value = int(10**-sigfigs)
    
    return sf_value, sf_error

#Data request
#---------------------------------------------------#
def mastQuery(request):
    global urlencode, httplib
    #"""Perform a MAST query.
    #    Parameters
    #    ----------
    #    request (dictionary): The MAST request json object
    #    Returns head,content where head is the response HTTP headers, and content is the returned data"""
    
    server='mast.stsci.edu'

    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    #Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)
    
    #opening the https connection
    conn = httplib.HTTPSConnection(server)

    #Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    #Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    # Close the https connection
    conn.close()
    
    return head, content

def catalogsearch_tic():
  global outData_tic, objRa, objDec
  mastRequest_tic = {"service":"Mast.Catalogs.Tic.Cone",
                 'params':{'ra':objRa,            #uses the planet ra and dec to search the tic catalog
                           'dec':objDec,
                           'radius':0.2},
                 'format':'json',
                 "timeout":10} 
  
  
  headers,outString = mastQuery(mastRequest_tic)
  outData_tic = json.loads(outString)
  
  print(outData_tic.keys())
  print("Query status:",outData_tic['status'])

def catalogsearch_gaia():
  global outData_gaia
  mastRequest_gaia = {'service':'Mast.Catalogs.GaiaDR2.Cone',
               'params':{'ra':objRa,              #uses planet ra and dec to search gaia catalog
                         'dec':objDec,
                         'radius':0.2},
               'format':'json',
               'page':5}
  
  headers,outString = mastQuery(mastRequest_gaia)
  outData_gaia = json.loads(outString)

  print(outData_gaia.keys())
  print("Query status:",outData_gaia['status'])

#Parameter calculations
#---------------------------------------------------#
#THIS IS THE MASS MULTIPLED BY sin(i)
def calc_mass():
  global p_mass, p_mass_err1, p_mass_err2, p_class, rounded_pmass, rounded_pmass_err1, rounded_pmass_err2
  X = ((period)/(2*np.pi*G))**(1/3)
  Y = 1-(ecc**2)
  Z = v_star*(M**(2/3))*np.sqrt(Y)  
  p_mass  = (X*Z)
  ecc_err2 = abs(neg_ecc_err2)        #PLANETARY MASS CALCULATION                                                                                             
  try:
    ecc_div_ecc1 = ecc_err1 / ecc
    ecc_div_ecc2 = ecc_err2 / ecc
  except:
    ecc_div_ecc1 = 0
    ecc_div_ecc2 = 0
  p_mass_err1 = p_mass*np.sqrt((v_star_err/v_star)**(2)+(period_err/(3*period))**(2)+((ecc_div_ecc1)**(2))+((2*M_err)/(3*M))**(2))    #PLANETARY MASS ERROR upper CALCULATION 
  p_mass_err2 = p_mass*np.sqrt((v_star_err/v_star)**(2)+(period_err/(3*period))**(2)+((ecc_div_ecc2)**(2))+((2*M_err)/(3*M))**(2))    #PLANETARY MASS ERROR lower CALCULATION 
  pmass = p_mass/J_mass             #CONVERSION OF MASS TO JUPITER MASSES
  pmass_err1 = p_mass_err1/J_mass
  pmass_err2 = p_mass_err2/J_mass     #CONVERSION OF MASS ERROR TO JUPITER MASSES
  rounded_pmass, errors = significant_figures(pmass, [pmass_err1,pmass_err2])
  rounded_pmass_err1 = errors[0]
  rounded_pmass_err2 = errors[1]             #ROUNDS THE MASS ERROR TO 1 SIG FIG
  print("Planetary mass =", rounded_pmass, "with errors of +", rounded_pmass_err1, "and -", rounded_pmass_err2 ,"Jupiter masses")      #PRINTS MASS
  classification = {1: 'Asteroidan',
         2: 'Mercurian',
         3: 'Subterran',
         4: 'Terran',
         5: 'Superterran',
         6: 'Neptunian',
         7: 'Jovian',
         8: 'bigger than Jovian'
         }
  pmass_earths = p_mass/Earth_mass
  #Checks what classification the planet is 
  if 0 <= pmass_earths < 0.00001:
    print("The planet is", classification[1])
    p_class = classification[1]
  elif 0.00001 <= pmass_earths < 0.1:
    print("The Planet is", classification[2])
    p_class = classification[2]
  elif 0.1 <= pmass_earths < 0.5:
    print("The Planet is", classification[3])
    p_class = classification[3]
  elif 0.5 <= pmass_earths < 2:
    print("The Planet is", classification[4])
    p_class = classification[4]
  elif 2 <= pmass_earths < 10:
    print("The Planet is", classification[5])
    p_class = classification[5]
  elif 10 <= pmass_earths < 50:
    print("The Planet is", classification[6])
    p_class = classification[6]
  elif 50 <= pmass_earths < 5000:
    print("The Planet is", classification[7])
    p_class = classification[7]
  elif pmass_earths >= 5000:
    print("The Planet is", classification[8])
    p_class = classification[8]


def calc_density():
    global rounded_pdens, rounded_pdens_err1, rounded_pdens_err2
    if p_radius ==0:
        rounded_pdens =0
        rounded_pdens_err1 =0
        rounded_pdens_err2 =0
        print("No density value")
    else:
      p_density = (p_mass/((4/3)*np.pi*(p_radius)**(3)))
      p_density_err1 = p_density*np.sqrt(((3*p_radius_err)/p_radius)**(2)+(p_mass_err1/(p_mass))**(2))
      p_density_err2 = p_density*np.sqrt(((3*p_radius_err)/p_radius)**(2)+(p_mass_err2/(p_mass))**(2))            #p density error   
      pdens = p_density/1000
      pdens_err1 = p_density_err1/1000 
      pdens_err2 = p_density_err2/1000    
      print(pdens_err1,pdens_err2)
      rounded_pdens, errors = significant_figures(pdens, [pdens_err1,pdens_err2])
      rounded_pdens_err1 = errors[0]    
      rounded_pdens_err2 = errors[1]
      print("Planetary density =", rounded_pdens, "with errors of +", rounded_pdens_err1, "and -", rounded_pdens_err2 , "g/cm3")

def calc_orbitalradius():
  global o_radius, o_radius_err, orad, orad_err, rounded_orad, rounded_orad_err
  o_radius = (((((period*24*3600)**2)*G*M)/(4*(np.pi)**2))**(1/3))
  o_radius_err = o_radius*(1/3)*np.sqrt(((2*period_err)/period)**(2)+(M_err/M)**(2))
  orad = o_radius/AU
  orad_err = o_radius_err/AU
  rounded_orad, rounded_orad_err = significant_figures(orad, orad_err)
  print("Orbital radius =", rounded_orad, "+/-", rounded_orad_err, "AU")

def calc_transittime():
  global rounded_dt, rounded_dt_err
  dt = (r_star*period*24)/(np.pi*o_radius)
  dt_err = dt*np.sqrt((period_err/period)**(2)+(r_star_err/r_star)**(2)+(o_radius_err/o_radius)**(2))
  rounded_dt, rounded_dt_err = significant_figures(dt, dt_err)
  print("Transit time =", rounded_dt, "+/-", rounded_dt_err, "hrs")

def calc_velocity():
  global rounded_pvel, rounded_pvel_err
  p_velocity = (2*np.pi*o_radius)/(period*24*3600*1000)
  p_velocity_err = p_velocity*np.sqrt((o_radius_err/o_radius)**(2)+(period_err/period)**(2))
  rounded_pvel, rounded_pvel_err = significant_figures(p_velocity, p_velocity_err)
  print("Planetary velocity =",rounded_pvel ,"+/-", rounded_pvel_err, "km/s")

def calc_pradius():
  global p_radius, p_radius_err, rounded_prad, rounded_prad_err
  p_radius = (r_star*np.sqrt(Drop_in_flux))
  p_radius_err = p_radius*np.sqrt((r_star_err/r_star)**(2)+(Drop_in_flux_error/(2*Drop_in_flux))**(2))
  prad = p_radius/J_radius
  prad_err = p_radius_err/J_radius
  try:
      rounded_prad, rounded_prad_err = significant_figures(prad, prad_err)
  except:
      rounded_prad = 0
      rounded_prad_err = 0
  print("Planetary radius =", rounded_prad, "+/-", rounded_prad_err, "Jupiter radii")

def calc_period_and_drop():
  global rounded_period, rounded_period_err, rounded_drop, rounded_drop_err
  rounded_period, rounded_period_err = significant_figures(period, period_err)
    
  rounded_drop, rounded_drop_err = significant_figures(Drop_in_flux, Drop_in_flux_error)
      
  print("Orbital period =", rounded_period, "+/-", rounded_period_err, "d")
  print("Drop in Flux =",rounded_drop,"+/-",rounded_drop_err,"normalised flux")
  return period

def calc_surfacetemp():
  global T_surf, T_surf_err, rounded_tsurf, rounded_tsurf_err, rounded_uhab, rounded_lhab, Hab_zone, lower_hab_err, upper_hab_err, rounded_uhab_err, rounded_lhab_err
  T_surf = (L/(16*np.pi*stef_bolt*(o_radius)**(2)))**(1/4)
  T_surf_err = T_surf*np.sqrt(((0.25*L_err/L)**(2))+((0.5*o_radius_err/o_radius)**(2)))
  rounded_tsurf, rounded_tsurf_err = significant_figures(T_surf, T_surf_err)
  #finds upper and lower limit of habitable zone
  upper_hab = ((273/T_surf)**(-2))*orad
  upper_hab_err = upper_hab*np.sqrt(((2*T_surf_err/T_surf)**(2))+((o_radius_err/o_radius)**(2)))
  lower_hab = ((373/T_surf)**(-2))*orad
  lower_hab_err = lower_hab*np.sqrt(((2*T_surf_err/T_surf)**(2))+((o_radius_err/o_radius)**(2)))
  rounded_uhab, rounded_uhab_err = significant_figures(upper_hab, upper_hab_err)
  rounded_lhab, rounded_lhab_err = significant_figures(lower_hab, lower_hab_err)

  print("The habitable zone of this star is between", rounded_lhab,"+/-",rounded_lhab_err, "and", rounded_uhab, "+/-", rounded_uhab_err, "AU")
  #checks whether the surface temperature is between 273K and 373K (freezing and boiling point of water respectively)
  if 273 <= T_surf <= 373:
   Hab_zone = 'Yes'
   print("Within habitable zone:", Hab_zone)
  else:
   Hab_zone = 'No'
   print("Within habitable zone:", Hab_zone)
   print("Surface temperature =", rounded_tsurf, "+/-", rounded_tsurf_err, "K")
   
def calc_grav():
  global rounded_g, rounded_g_err1, rounded_g_err2
  #Calculating the surface gravity of the planet
  if p_mass == 0:
      rounded_g =0
      rounded_g_err1 =0
      rounded_g_err2 =0
      print("No Surface gravity value")
  else:   
      try:
          g = ((G*p_mass)/(p_radius)**(2))
          g_err1 = g*np.sqrt(((p_mass_err1/p_mass)**(2))+(((2*p_radius_err)/p_radius)**(2)))
          g_err2 = g*np.sqrt(((p_mass_err2/p_mass)**(2))+(((2*p_radius_err)/p_radius)**(2)))
          rounded_g, errors = significant_figures(g, [g_err1,g_err2])
          rounded_g_err1 = errors[0]
          rounded_g_err2 = errors[1]
          print("Surface gravity =", rounded_g, "+/-", rounded_g_err1, rounded_g_err2, "m/s^2")
      except:
          rounded_g =0
          rounded_g_err1 =0
          rounded_g_err2 =0     

  if p_radius == 0:
      rounded_g =0
      rounded_g_err1 =0
      rounded_g_err2 =0
      print("No Surface gravity value")
  else:   
      try:
          g = ((G*p_mass)/(p_radius)**(2))
          g_err1 = g*np.sqrt(((p_mass_err1/p_mass)**(2))+(((2*p_radius_err)/p_radius)**(2)))
          g_err2 = g*np.sqrt(((p_mass_err2/p_mass)**(2))+(((2*p_radius_err)/p_radius)**(2)))
          rounded_g, errors = significant_figures(g, [g_err1,g_err2])
          rounded_g_err1 = errors[0]
          rounded_g_err2 = errors[1]
          print("Surface gravity =", rounded_g, "+/-", rounded_g_err1, rounded_g_err2, "m/s^2") 
      except:
          rounded_g =0
          rounded_g_err1 =0
          rounded_g_err2 =0

      
 
#Jordan's Main code
#---------------------------------------------------#
def QueryParameters(Star):
    global urlencode, httplib, objRa, objDec, outData_tic, outData_gaia, exoarchive
    
    try: # Python 3.x
        from urllib.parse import quote as urlencode
        from urllib.request import urlretrieve
    except ImportError:  # Python 2.x
        from urllib import pathname2url as urlencode
        from urllib import urlretrieve

    try: # Python 3.x
        import http.client as httplib 
    except ImportError:  # Python 2.x
        import httplib
    
    pp = pprint.PrettyPrinter(indent=4)
    resolverRequest = {'service':'Mast.Name.Lookup','params':{'input':Star,'format':'json'},}
    
    headers,resolvedObjectString = mastQuery(resolverRequest)
    resolvedObject = json.loads(resolvedObjectString)

    #pp.pprint(resolvedObject)

    #Extracts the ra and dec of planet
    objRa = resolvedObject['resolvedCoordinate'][0]['ra']
    objDec = resolvedObject['resolvedCoordinate'][0]['decl']

    cachebreaker = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    
    #searches MAST for 2 different catalogs
    
    catalogsearch_gaia()
    catalogsearch_tic()
    
    #Stellar RA and Dec needed to search the exoplanet archive (NASA)
    stellarRA = round(outData_tic['data'][0]['RA_orig'], 7)
    stellarDEC = round(outData_tic['data'][0]['Dec_orig'], 7)
    
    
    #EXOPLANET ARCHIVE QUERY to find the eccentricity values of the planet needed to calculate mass
    service = vo.dal.TAPService("https://exoplanetarchive.ipac.caltech.edu/TAP")
    SQL = "SELECT pl_orbeccen, pl_orbeccenerr1, pl_orbeccenerr2, rowupdate FROM ps WHERE ra = " + str(stellarRA) + " AND dec = " + str(stellarDEC) + "AND pl_orbeccen >= 0 AND pl_orbeccenerr1 IS NOT NULL ORDER BY rowupdate DESC "
    exoarchive = service.search(SQL)
    
def CalculateParameters(Star):
    global rounded_pmass,rounded_pmass_err1,rounded_pmass_err2,rounded_pdens,rounded_pdens_err1,rounded_pdens_err2,rounded_tsurf,rounded_tsurf_err, rounded_g, rounded_g_err1, rounded_g_err2
    global rounded_orad,rounded_orad_err,rounded_dt,rounded_dt_err,rounded_pvel,rounded_pvel_err,rounded_prad,rounded_prad_err,rounded_period,rounded_period_err,p_radius,p_radius_err,p_mass,p_mass_err1,p_mass_err2
    global classification, p_class, Hab_zone, second_attempt, star_variables, M, L, r_star, v_star, M_err, L_err, v_star_err, r_star_err
    #Constants needed to convert values to Astronomical units
    global J_mass,J_radius,AU,Solar_radius,Solar_mass,G,solar_L,stef_bolt,Earth_mass
    J_mass = 1.898e27
    J_radius = 69911000
    AU = 1.5e11
    Solar_radius = 6.96e8
    Solar_mass = 1.989e30
    G = 6.67e-11
    solar_L = 3.83e26
    stef_bolt = 5.67e-8
    Earth_mass = 5.972e24
    
    #PARAMETERS USED FOR ROUNDING
    global sig_fig_1, sig_fig_2
    sig_fig_1 = 3
    sig_fig_2 = 1
    
    #PARAMETERS FOUND USING THE QUERIES AND CATALOGS
    global r_star,r_star_err,v_star,v_star_err,M,M_err,ecc,ecc_err1,neg_ecc_err2,L,L_err
    if second_attempt == False:
        r_star = outData_tic['data'][0]['rad'] #in m
        v_star = outData_gaia['data'][0]['radial_velocity']  #in km/s
        M = outData_tic['data'][0]['mass']  #in kg
        L = outData_tic['data'][0]['lum']
        M_err = outData_tic['data'][0]['e_mass']  #in kg
        v_star_err = outData_gaia['data'][0]['radial_velocity_error']  #in km/s
        r_star_err = outData_tic['data'][0]['e_rad']   #in m
        L_err = outData_tic['data'][0]['e_lum']
        star_variables = [M,r_star,v_star,L]
        
        if M_err is None:
          M_err =0
        else:
          M_err = M_err*Solar_mass
        if v_star_err is None:
          v_star_err =0
        else:
          v_star_err = abs(v_star_err)*1000
        if L_err is None:
          L_err =0
        else:
          L_err = L_err*solar_L
        if r_star_err is None:
          r_star_err =0
        else:
          r_star_err = r_star_err*Solar_radius
    
    if len(exoarchive) == 0:
      ecc = None
      ecc_err1 =None
      neg_ecc_err2 = None
    else:
      ecc = exoarchive[0]['pl_orbeccen']
      ecc_err1 = exoarchive[0]['pl_orbeccenerr1']
      neg_ecc_err2 = exoarchive[0]['pl_orbeccenerr2']
    
    if math.isnan(period) and math.isnan(period_err) or period == 0 and period_err ==0:
        print("No Values")
        #sets all values to zero for the data frame
        NoResults(Star)
        return

    else:
        if ecc is None and ecc_err1 is None and neg_ecc_err2 is None:
          ecc=0
          ecc_err1 = 0
          neg_ecc_err2 =0
          calc_period_and_drop()
          if M is None:
            print("No Values")
            #sets all values to zero
            NoResults(Star)
            return
          else:
            M = M*Solar_mass
            calc_orbitalradius()
            calc_velocity()
            if v_star is None:
              print("No Planetary mass")
              p_mass =0
              p_mass_err1 =0
              p_mass_err2 =0
              rounded_pmass = 0
              rounded_pmass_err1 = 0
              rounded_pmass_err2 =0
              p_class = None
            else:
              v_star = abs(v_star)*1000
              calc_mass()
            if L is None:
                print("No Surface Temperature")
                rounded_tsurf =0
                rounded_tsurf_err = 0
            else:
              L = L*solar_L
              calc_surfacetemp()
            if math.isnan(Drop_in_flux) and math.isnan(Drop_in_flux_error) or Drop_in_flux ==0 and Drop_in_flux_error ==0:
              print("No Planetary radius, No Planetary density")
              rounded_prad= 0
              rounded_prad_err= 0
              rounded_pdens = 0
              rounded_pdens_err1= 0
              rounded_pdens_err2 =0
              
            else:
              if r_star is None:
                print("No transit time, No Planetary radius, No Planetary density, No surface gravity")
                rounded_prad= 0
                rounded_prad_err= 0
                rounded_pdens = 0
                rounded_pdens_err1 =0
                rounded_pdens_err2 =0
                rounded_dt =0
                rounded_dt_err =0
                rounded_g_err1 =0
                rounded_g_err2 =0
                rounded_g =0
              else:
                r_star = r_star*Solar_radius
                calc_transittime()
                calc_pradius()
                calc_grav()
                if v_star is None:
                  print("No Planetary density")
                  rounded_pdens = 0
                  rounded_pdens_err1 =0
                  rounded_pdens_err2 =0
                  
                else:
                  calc_density()
        else:
          calc_period_and_drop()
          if M is None:
            print("No Values")
            #sets all values to zero
            NoResults(Star)
            return
          else:
            M = M*Solar_mass
            calc_orbitalradius()
            calc_velocity()
            if v_star is None:
              print("No Planetary mass")
              p_mass =0
              p_mass_err1 =0
              p_mass_err2 =0

              rounded_pmass = 0
              rounded_pmass_err1 =0
              rounded_pmass_err2 =0
              p_class = None
            else:
              v_star = abs(v_star)*1000
              calc_mass()
            if L is None:
              print("No Surface Temperature")
              rounded_tsurf =0
              rounded_tsurf_err = 0
            else:
              L = L*solar_L
              calc_surfacetemp()
            if math.isnan(Drop_in_flux) and math.isnan(Drop_in_flux_error) or Drop_in_flux ==0 and Drop_in_flux_error ==0:
              print("No Planetary radius, No Planetary density")
              rounded_prad= 0
              rounded_prad_err= 0
              rounded_pdens = 0
              rounded_pdens_err1 = 0
              rounded_pdens_err2 =0
            else:
              if r_star is None:
                print("No transit time, No Planetary radius, No Planetary density, No surface gravity")
                rounded_prad= 0
                rounded_prad_err= 0
                rounded_pdens = 0
                rounded_pdens_err1=0
                rounded_pdens_err2 =0
                rounded_dt =0
                rounded_dt_err =0
                rounded_g= 0
                rounded_g_err1 =0
                rounded_g_err2 =0
              else:
                r_star = r_star*Solar_radius
                calc_transittime()
                calc_pradius()
                calc_grav()
                if v_star is None:
                  print("No Planetary density")
                  rounded_pdens = 0
                  rounded_pdens_err1 =0
                  rounded_pdens_err2 =0
                else:
                  calc_density()
            
    SaveResults(Star)
    
def NoResults(Star):
    global rounded_pmass,rounded_pmass_err1,rounded_pmass_err2,rounded_pdens,rounded_pdens_err1,rounded_pdens_err2, rounded_tsurf,rounded_tsurf_err,rounded_g_err1 ,rounded_g_err2, rounded_g
    global rounded_orad,rounded_orad_err,rounded_dt,rounded_dt_err,rounded_pvel,rounded_pvel_err,rounded_prad,rounded_prad_err,rounded_period,rounded_period_err, rounded_drop, rounded_drop_err
    global p_class,Hab_zone, period, period_err, Drop_in_flux, Drop_in_flux_error
    if period != 0:
        rounded_period, rounded_period_err = significant_figures(period, period_err)
        rounded_drop, rounded_drop_err = significant_figures(Drop_in_flux, Drop_in_flux_error)
    else:
        rounded_period= 0
        rounded_period_err= 0
        rounded_drop= 0
        rounded_drop_err= 0
        
    print("No available values") 
    p_class = None
    rounded_pmass = 0
    rounded_pmass_err1 = 0
    rounded_pmass_err2 = 0
    rounded_pdens = 0
    rounded_pdens_err1 = 0
    rounded_pdens_err2 = 0
    rounded_orad= 0
    rounded_orad_err= 0
    rounded_dt= 0
    rounded_dt_err= 0
    rounded_pvel= 0
    rounded_pvel_err= 0
    rounded_prad= 0
    rounded_prad_err= 0
    rounded_tsurf= 0
    rounded_tsurf_err= 0
    rounded_g =0
    rounded_g_err1 =0
    rounded_g_err2 =0
    Hab_zone='No'
    
    SaveResults(Star)

#-------------------------------------------------------------------------#
#Will GUI code 
#-------------------------------------------------------------------------# 

def User_values(master,text,value):
    frame = LabelFrame(master,borderwidth=0)
    title = Label(frame,text=text[0],width=10)
    title.pack(side='left')
    title.pack_propagate(0)
    input_box = Entry(frame, bd=8,font='{Arial Nova Cond Light} 14',width=30)
    units = Label(frame,text=text[1])
    if value != None:
        input_box.insert('end', value)
    input_box.pack(side='left',padx=10)
    units.pack(side='left')
    frame.pack(anchor='nw',padx=20,pady=5)
    
def CheckForMissingValues(star_titles,star_variables):
    global T_surf, in_tsurf, waiting
    input_values = messagebox.askquestion('Warning','Some of the values needed to calculate the exoplanet parameters could not be found. \
Do you want to input values yourself? Answering NO will cause the program to return the values it \
currently has.')
    if input_values == 'yes':
        user_win = Toplevel()
        user_win.wm_title("Input Values")
        user_win.geometry("900x550")
        
        screen_width = user_win.winfo_screenwidth()
        screen_height = user_win.winfo_screenheight()
        user_win.geometry("+"+str(int(screen_width/2)-450)+"+"+str(int(screen_height/2)-275))
        
        title = Label(user_win,text='STAR VARIABLES',fg='darkblue', font='{Arial Nova Cond Light} 25 underline')
        title.pack(anchor='center',padx=20,pady=(10,5))
        How_To = Label(user_win,text='Input values for any blank fields. You do not need to fill all \
of them but having values for each will yield the maximum amount of results. The star Mass is the most import value \
as no calculations can be done without it.',justify='left',wraplengt=860)
        How_To.pack(anchor='nw',padx=20,pady=(0,10))
        
        for i in range(len(star_variables)):
            User_values(user_win,star_titles[i],star_variables[i])
        waiting = IntVar()
        cancel_btn = Button(user_win, text='Cancel',command=lambda: waiting.set(2))
        cancel_btn.pack(side='right',padx=(5,20),pady=20,anchor='ne')
        submit_btn = Button(user_win, text='Submit',command=lambda: letsgo(user_win))
        submit_btn.pack(side='right',pady=20,anchor='ne')
        user_win.protocol("WM_DELETE_WINDOW",lambda: waiting.set(2))
        
        submit_btn.wait_variable(waiting)               
        user_win.destroy()
        
        if waiting == 2:
            return False
        else:
            return True
        
    elif input_values == 'no':
        return False
   
def letsgo(user_win):
    global waiting, M, r_star, v_star,  L, M_err, L_err, v_star_err, r_star_err
    values = [child.winfo_children()[1].get() for child in user_win.winfo_children() if child.winfo_class() == 'Labelframe']
    numeric_values = [CheckIfNumber(val) for val in values]
   
    if 'Non-Value' in numeric_values:
        messagebox.showwarning('Warning','Values must be numeric!')
        return
    else:
        M = numeric_values[0]
        r_star = numeric_values[1]
        v_star = numeric_values[2]
        L = numeric_values[3]
        if M != None: M_err = M*0.05*Solar_mass
        if v_star != None: v_star_err = v_star*0.05*1000
        if r_star != None: r_star_err = r_star*0.05*Solar_radius
        if L != None: L_err = L*0.05*solar_L
        waiting.set(1)
        return
    
def CheckIfNumber(val):
    if val.isnumeric() == True:
        return int(val)
    if val == '':
        return None
    else:
        try:
            float(val)
            return float(val)
        except:
            return 'Non-Value'
#-------------------------------------------------------------------------#
#Parin's code 
#-------------------------------------------------------------------------#    

def SaveResults(Star):
    global data_frame, mode, Final_Results, M, r_star, v_star,  L, M_err, L_err, v_star_err, r_star_err, second_attempt, star_variables
    star_titles = [['Mass:','M☉ (Solar Mass)'],['Radius:','R☉ (Solar Radius)'],['Velocity:','km/s'],['Luminosity:','L☉ (Solar Luminosity)']]
    
    if mode == 1 and None in star_variables and second_attempt == False and rounded_period != 0:
        second_attempt = True
        rerun = CheckForMissingValues(star_titles,star_variables)
        if rerun == True and (M != None or r_star != None or v_star != None or L != None):
            CalculateParameters(Star)
            return
            
    data_frame = data_frame.append({'Star': Star,
                                    'Star Radius': r_star,
                                    'Planet Classification': p_class,
                                    'Planetary mass (M⊕)': rounded_pmass,
                                    'Mass Error1': rounded_pmass_err1,
                                    'Mass Error2': rounded_pmass_err2,
                                    'Planetary density (g/cm\u00b3)': rounded_pdens, 
                                    'Density Error1': rounded_pdens_err1,
                                    'Density Error2': rounded_pdens_err2,
                                    'Orbital radius (AU)': rounded_orad,
                                    'Orbital radius Error': rounded_orad_err,
                                    'Transit time (hrs)': rounded_dt, 
                                    'Transit Error': rounded_dt_err,
                                    'Planetary velocity (km/s)': rounded_pvel,
                                    'Velocity Error': rounded_pvel_err,
                                    'Planetary radius (R⊕)': rounded_prad, 
                                    'Planetary radius Error': rounded_prad_err,
                                    'Drop in flux': rounded_drop,
                                    'Drop in flux Error': rounded_drop_err,
                                    'Period (d)': rounded_period,
                                    'Period Error': rounded_period_err,
                                    'Surface temperature (K)': rounded_tsurf,
                                    'Surface temperature Error': rounded_tsurf_err,
                                    'Surface gravity (m/s\u00b2)': rounded_g,
                                    'Surface gravity Error1': rounded_g_err1,
                                    'Surface gravity Error2': rounded_g_err2,
                                    'Within habitable zone': Hab_zone}, ignore_index=True)
    
    if M is not None: M = M/Solar_mass
    if L is not None: L = L/solar_L
    if r_star is not None: r_star = r_star/Solar_radius
    if v_star is not None : v_star = v_star/1000
    
    if second_attempt == False:
        if M_err is not None: M_err = M_err/Solar_mass
        if L_err is not None: L_err = L_err/solar_L
        if r_star_err is not None: r_star_err = r_star_err/Solar_radius
        if v_star_err is not None : v_star_err = v_star_err/1000
        
    return



