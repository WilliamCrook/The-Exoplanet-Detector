"""CODE AUTHOR: PARIN VYAS"""
from vpython import vector, sphere, mag, local_light, color, vec, ring, label, gcurve, rate, scene, graph
import numpy as np


def MainProgram(Final_Results):
    exoplanet_variables = []  
    
    ##  Constants
    # Radii:
    starrad = Final_Results['Star Radius'].iloc[0]
    mercuryrad = 2.44e6
    venusrad = 6.05e6
    earthrad = 6.38e6
    marsrad = 3.4e6
    
    for ind in range(len(Final_Results)):
        exoplanet_variables.append([
                            Final_Results['Planetary radius (R⊕)'].iloc[ind]*earthrad,
                            Final_Results['Orbital radius (AU)'].iloc[ind]*1.496e11,
                            Final_Results['Period (d)'].iloc[ind]*24,
                            ])
    
    # Position:
    mercurypos = vector(0, 0, (5.79e10+starrad))
    venuspos = vector(0, 0, (1.082e11+starrad))
    earthpos = vector(0, 0, (1.496e11+starrad))
    marspos = vector(0, 0, (2.279e11+starrad))
    
    # Orbital Periods (Hours):#
    mercuryperiod = (88*24)
    venusperiod = (224.7*24)
    earthperiod = (365.2*24)
    marsperiod = (687*24)
    
    # CHZ Attributes
    chz_inner = ((0.85*mag(earthpos)) + starrad)
    chz_outer = ((1.6*mag(earthpos)) + starrad)
    chz_thickness = (0.009*mag(earthpos))
    
    # Time
    t = 0
    deltat = 0.1
    
    # Scene and Lights
    star_light = local_light(pos=vector(0,0,0), color=color.white, visible = True)
    
    
        ##  Objects
    star=sphere(color=color.orange,pos=vec(0,0,0),radius=starrad, emissive=True, texture = "http://i.imgur.com/yoEzbtg.jpg")
    mercury=sphere(color=color.gray(0.75),pos=mercurypos,radius=mercuryrad, shininess=10, make_trail=True, retain=5000, opacity = 0.3)
    venus=sphere(color=color.gray(0.75),pos=venuspos,radius=venusrad, shininess=10, make_trail=True, retain=5000,opacity = 0.3)
    earth=sphere(color=color.gray(0.75),pos=earthpos,radius=earthrad, shininess=10, make_trail=True,retain=5000, opacity = 0.3)
    mars=sphere(color=color.gray(0.75),pos=marspos,radius=marsrad, shininess=10, make_trail=True,retain=5000, opacity = 0.3)
    ring1 = ring(pos= vec(0,0,0), axis = vec(0,1,0), radius = chz_inner, thickness = chz_thickness, color = color.green, opacity = 0.75)
    ring2 = ring(pos= vec(0,0,0), axis = vec(0,1,0), radius = chz_outer, thickness = chz_thickness, color = color.red, opacity = 0.75)
    
    exoplanets = []
    for exo in exoplanet_variables:
        exoplanet=sphere(color=color.cyan,pos=vector(0, 0,(exo[1] + starrad)),radius=exo[0], shininess=10, make_trail=True,retain=5000)
        exoplanet_label = label(pos=vector(0, 0, mag(exoplanet.pos)),
            text='Exoplanet', xoffset=10, color = color.cyan,
            line = False, yoffset=10, height=12, 
            border=3, font='sans')
        exoplanets.append([exoplanet,exoplanet_label])
    
        ##  Labels
    star_label = label(pos=star.pos,
        text='Star', xoffset=10, opacity = 0.5,
        line = False, yoffset=10, height=12, 
        border=3, font='sans')
    
    mercury_label = label(pos=vector(0, 0, mag(mercury.pos)),
        text='Mercury', xoffset=10, opacity = 0.5,
        line = False, yoffset=10, height=12, 
        border=3, font='sans')
    
    venus_label = label(pos=vector(0, 0, mag(venus.pos)),
        text='Venus', xoffset=10, opacity = 0.5,
        line = False, yoffset=10, height=12, 
        border=3, font='sans')
    
    earth_label = label(pos=vector(0, 0, mag(earth.pos)),
        text='Earth', xoffset=10, opacity = 0.5,
        line = False, yoffset=10, height=12, 
        border=3, font='sans')
    
    mars_label = label(pos=vector(0, 0, mag(mars.pos)),
        text='Mars', xoffset=10, opacity = 0.5,
        line = False, yoffset=10, height=12, 
        border=3, font='sans')
    
    ring_label = label(pos=(ring1.pos+vector(ring1.radius, 0, 0)),
        text='CHZ Inner', xoffset=4, color = color.green,
        line = False, yoffset=12, height=12, 
        border=3, font='sans')
    
    ring_label2 = label(pos=(ring2.pos+vector(ring2.radius, 0, 0)),
        text='CHZ Outer', xoffset=4, color = color.red,
        line = False, yoffset=12, height=12, 
        border=3, font='sans')
    
    
        ##  Graphs
    gd = graph(width=635, height=600, xtitle='Time', ytitle='Position (Exoplanet)',
               title='\n<b>Graph showing Position of the Exoplanet</b>',
               foreground=color.gray(0.5), background=color.white,
               xmax=5000, xmin=0, ymax=5e10, ymin=-5e10)
    g1 = gcurve(color=color.red, label = "Z-Position")
    g2 = gcurve(color=color.cyan, label = "X-Position") # a graphics curve
    
    scene.camera.follow(star)
    
    while True:
        rate(10)
        scene.camera.follow(star)
        t = t + deltat
        star.rotate(angle=((2*np.pi)/365), axis=vector(0,1,0), origin=vector(0,0,0))
        for jnd in range(len(exoplanets)):
            exoplanets[jnd][0].rotate(angle=((2*np.pi)/exoplanet_variables[jnd][2]), axis=vector(0,1,0), origin=vector(0,0,0))
        mercury.rotate(angle=((2*np.pi)/mercuryperiod), axis=vector(0,1,0), origin=vector(0,0,0))
        venus.rotate(angle=((2*np.pi)/venusperiod), axis=vector(0,1,0), origin=vector(0,0,0))
        earth.rotate(angle=((2*np.pi)/earthperiod), axis=vector(0,1,0), origin=vector(0,0,0))
        mars.rotate(angle=((2*np.pi)/marsperiod), axis=vector(0,1,0), origin=vector(0,0,0))
        g1.plot((t, exoplanets[0][0].pos.z))
        g2.plot((t, exoplanets[0][0].pos.x))