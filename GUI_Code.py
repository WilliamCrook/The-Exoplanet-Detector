"""CODE AUTHOR: WILLIAM CROOK"""
#---------------------------------------------------------------------------------------------------------------#
"""IMPORTS"""
#---------------------------------------------------------------------------------------------------------------#

from astropy.io import fits
import os
import webbrowser
import threading
import requests
from datetime import datetime

import lightkurve as lk
import pandas as pd
import csv

from tkinter import (
    Button, Canvas, Entry, font, Frame, IntVar, Label, LabelFrame, messagebox,   
    OptionMenu, PhotoImage, StringVar, Tk, ttk 
    )
from tkinter.filedialog import askdirectory, askopenfilename 

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

#Scripts written for project#
#-----------------------------------------------------------------------------#
import Analysis_Code
import Parameter_Calculations
import Make_Model
#-----------------------------------------------------------------------------#

import warnings
warnings.filterwarnings("ignore")

import ctypes
ctypes.windll.shcore.SetProcessDpiAwareness(1)

#---------------------------------------------------------------------------------------------------------------#
"""DATA PROCESSING"""
#---------------------------------------------------------------------------------------------------------------#

def TpfToLightCurve(tpf):
    try:
        #Converts tpf to lightcurve by making a threshold mask#
        target_mask = tpf.create_threshold_mask(threshold=5, reference_pixel='center')
        target_lc = tpf.to_lightcurve(aperture_mask=target_mask)
        n_target_pixels = target_mask.sum()
        
        #Creates a mask to eliminate light from background pixels#
        background_mask = ~tpf.create_threshold_mask(threshold=0.001, reference_pixel=None)
        n_background_pixels = background_mask.sum()
        background_lc_per_pixel = tpf.to_lightcurve(aperture_mask=background_mask) / n_background_pixels
        background_estimate_lc = background_lc_per_pixel * n_target_pixels
        lc = target_lc - background_estimate_lc.flux
        
        #Cleans and flattens the lightcurve#
        lc_clean = lc.remove_outliers(sigma=20)
        lc_flat = lc_clean.flatten()
    except:
        #Returns None in the case of an error#
        return
    
    return lc_flat

def saveData(search_result,iterator,path,DataType):  
    global mode, lc_list    
    for i in range(len(iterator)):
        #Sets lc to None to stop an error being thrown#
        lc = None
        #Finds the sector number for the data#
        Sector = search_result[iterator[i]].table['sequence_number'][0]
        #Converts data to from tpf to light curve for FFI or downloads light curve for SPOC or TESS-SPOC#
        try:
            if DataType == 'FFI':
                tpf = search_result[i].download(cutout_size=10)
                lc = TpfToLightCurve(tpf)
            if DataType == 'SPOC' or DataType == 'TESS-SPOC':
                lc = search_result[iterator[i]].download()    
        except Exception as e:
            #Saves an txt file in case an error is thrown for future reference#
            file = open(str(path)+'/'+str(Sector)+'.txt','w+')
            file.write('Something went wrong.')
            file.write(str(e))
            file.close()
            
        if lc != None:
            #Converts file to fits and saves it in folder#
            exptime = search_result[iterator[i]].exptime.value[0]
            lc.to_fits(str(path)+'/'+str(Sector)+'.fits', overwrite=True)
            hdu = fits.open(str(path)+'/'+str(Sector)+'.fits',mode='update')
            hdu[0].header['exptime'] = exptime
            hdu.flush()
            lc_list.append([Sector,exptime,lc])
        else:
            #Creates a txt file in the case where no data was found for the sector#
            file = open(str(path)+'/'+str(Sector)+'.txt','w+')
            file.write('Due to an error downloading the data a light curve is not available for this sector.')
            file.close()
            
def CheckforNewData(Star,DataType,processing):
    global database_path
    new_download_check = 0
    lk.log.setLevel('INFO')
    #Searches for new data based on data type and if it finds some saves it in the correct folder#
    folderpath = database_path+'/'+str(Star)
    processing.config(text='Checking internet for new data...')
    if DataType == 'FFI' or DataType == 'ALL':
        search_result = lk.search_tesscut(Star)
        search_result = checkifexists(Star,search_result,folderpath+'/FFI')
        if len(search_result) != 0:
            new_download_check = 1
            processing.config(text='New data found for '+Star+' downloading...')
            findFolder(search_result,folderpath,"FFI")
        
    if DataType == 'SPOC' or DataType == 'ALL':
        search_result = lk.search_lightcurve(Star, author='SPOC',exptime=120)
        search_result = checkifexists(Star,search_result,folderpath+'/SPOC')
        if len(search_result) != 0:
            new_download_check = 1
            processing.config(text='New data found for '+Star+' downloading...')
            findFolder(search_result,folderpath,"SPOC")
            
    if DataType == 'TESS-SPOC' or DataType == 'ALL':
        search_result = lk.search_lightcurve(Star, author='TESS-SPOC')
        search_result = checkifexists(Star,search_result,folderpath+'/TESS-SPOC')
        if len(search_result) != 0:
            new_download_check = 1
            processing.config(text='New data found for '+Star+' downloading...')
            findFolder(search_result,folderpath,"TESS-SPOC")
            
    if new_download_check == 0 and DataType == 'ALL':
        processing.config(text='You already have all the data for this star.')
    elif new_download_check == 0:
        #Tells user when no data or a certain data type is available#
        print('TESS has no',DataType,'data available for this Star.')

def checkifexists(Star,search_result,folderpath):
    try:
        already_downloaded = [int(file.replace('.fits','')) for file in os.listdir(folderpath) if file.endswith('.fits')==True]
        newsectors = search_result.table['sequence_number'].tolist()
        removelist = [newsectors.index(x) for x in newsectors if x in already_downloaded]
        removelist.reverse()
        for i in range(len(removelist)): 
            search_result.table.remove_row(removelist[i])
        return search_result
    except:
        return search_result

def RunProgram(Star,currenttype,processing):
    global lc_list, dfgraph, runningpage, DataType, tv1, database_path, checked_stars, Results_list, MT_df, Final_Results
    lc_list = []
    dfgraph = []
    folderpath = database_path+'/'+str(Star)
    GetLightcurves(Star,currenttype,folderpath,processing)
    print(Star)
    if len(lc_list) == 0:
        if runningpage == 1:
            messagebox.showwarning('Warning','There is no data available for this star.')
        if runningpage == 2:
            MT_df = MT_df.append(pd.DataFrame([{'Star':Star,'Period (d)':'There is no data available for this star.'}]))
        return False
    
    try:
        dfgraph, Results = Analysis_Code.MainProgram(Star,lc_list,DataType,runningpage)
        Final_Results = Parameter_Calculations.MainProgram(Star,checked_stars,Results,runningpage)
    except Exception as e:
        print(e)
        if runningpage == 1:
            messagebox.showwarning('Warning','Something went wrong hun :(  -  '+str(e))
        if runningpage == 2:
            MT_df = MT_df.append(pd.DataFrame([{'Star':Star,'Period (d)':'An error occured :('}]))
        return False
    
    if len(Final_Results) != 0 and runningpage == 1:
        if Final_Results['Period (d)'].iloc[0] == 0:
            messagebox.showerror('Warning','Did not detect any exoplanets for star.')
            return False
        Final_Results = Final_Results.T
        Results_list = []
        for col in Final_Results:
            formatted_results = pd.DataFrame(columns=['Parameter','Value','Error'])
            
            df = Final_Results[[col]]
            df = df.rename(columns={col: 0})
            
            for i in range(len(df)):
                header = df.iloc[[i]].index[0]
                value = df.iloc[[i]][0][0]
                if header == 'Star' or header == 'Star Radius' or header.endswith('Error') == True or header.endswith('Error',0,-1) == True:
                    continue
                if i == len(df)-1:
                    error = ''
                else:
                    error_header = df.iloc[[i+1]].index[0]
                    if error_header.endswith('Error') == True:
                        error = '±',df.iloc[[i+1]][0][0]
                    elif error_header.endswith('Error',0,-1) == True:
                        error = '+',df.iloc[[i+1]][0][0],'-',df.iloc[[i+2]][0][0]
                    else:
                        error = ''
                
                formatted_results = formatted_results.append({'Parameter':header,
                                                              'Value':value,
                                                              'Error':error
                                                            }, ignore_index=True)
            Results_list.append(formatted_results)
            
        MakeTable(Results_list[0]) 
        Final_Results = Final_Results.T

        
    elif runningpage == 1:
        messagebox.showerror('Warning','Something went wrong.')
        
    elif runningpage == 2:
        MT_df = MT_df.append(Final_Results)
        
    return True

#---------------------------------------------------------------------------------------------------------------#
"""DATABASE SEARCHING"""
#---------------------------------------------------------------------------------------------------------------#

def findFolder(search_result,folderpath, DataType):
    #Checks if a folder for this star already exists and if not creates one#
    if os.path.exists(folderpath) == False:
        try:
            os.mkdir(folderpath)
            os.mkdir(folderpath+'/FFI')
            os.mkdir(folderpath+'/SPOC')
            os.mkdir(folderpath+'/TESS-SPOC')
        except OSError:
            print('Could not create save folder.')
        else:
            #Save data to the folder#
            saveData(search_result,range(len(search_result)),folderpath+'/'+str(DataType),DataType)
    else:
        #If the folder exists it saves the data in it#
        saveData(search_result,range(len(search_result)),folderpath+'/'+str(DataType),DataType)
        
def GetLightcurves(Star,currenttype,folderpath,processing):
    global filetype, DataType
    #Selects the datatype in a format that easy to read#
    DataType = filetype[currenttype]
    
    if os.path.exists(folderpath) == False:
        #Searches for new data for a valid star name when the folder does not already exist#
        print('No Data in database checking internet')
        CheckforNewData(Star,DataType,processing)
    else:
        #Finds the folder of the data type required and call the light curves#
        path = folderpath+'/'+str(DataType)
        #Makes sure to only return fits file type rather than the txt files of failed sectors#
        Sectors = [file for file in os.listdir(str(path)) if file.endswith(".fits") == True]
        CallLightcurves(Sectors,path,Star,DataType,currenttype,processing)
            
def CallLightcurves(Sectors,path,Star,DataType,currenttype,processing):
    global filetype, lc_list
    #Recovers the sectors from the file names#
    sector_list = sorted([int(os.path.splitext(Sectors[i])[0]) for i in range(len(Sectors))])
        
    if len(sector_list) == 0:
        print("There is no",filetype[currenttype],"data available for this star.")
        print("Checking internet for new",DataType,'data.')
        CheckforNewData(Star,DataType,processing)
    
    for i in range(len(sector_list)):
        #Returns all the fits files in the folder by sector#
        Sector = sector_list[i]
        filename = str(path)+'/'+str(Sector)+'.fits'
        hdu = fits.open(filename,mode="readonly")
        
        with hdu as hdulist:      
            exptime = hdulist[0].header['exptime']
            times = hdulist[1].data['time']
            pdcsap_fluxes = hdulist[1].data['FLUX']
            errors = hdulist[1].data['FLUX_ERR']
        
        #Converts the hdu file into light curves#
        lc = lk.LightCurve(time=times,flux=pdcsap_fluxes,flux_err=errors)
        lc_list.append([Sector,exptime,lc])  

def FindName(userinput):
    global dictionary_df, misc_path
    
    formattedinput = userinput.lower().replace("-","").replace(" ","")
        
    #Searches the dictionary for a match and returns row number of match#
    searchdict = list(dictionary_df.isin([str(formattedinput)]).any(1).values.tolist())
    row = [searchdict.index(x) for x in searchdict if x == True]
    if len(row) == 0:
        tic_check = []
        response_df = pd.DataFrame()
        userinput = userinput.upper()
        first_result = None
        
        try:
            base_url = "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=aliastable&objname="
            response = requests.get(base_url+userinput+'&format=json').json()
            response_df = response_df.append(response,ignore_index=True)
            tic_check = [x for x in response_df['aliasdis'] if str(x).startswith('TIC')]
            first_result = response_df.iloc[0][0]
        except:
            pass
        
        if len(tic_check) != 0 and userinput.startswith('TIC') != True:
            temp_df = pd.DataFrame([{'aliasdis':tic_check[0]}])    
        elif userinput.startswith('TIC'):
            squished_Star = userinput.replace(" ","")
            temp_df = pd.DataFrame([{'aliasdis':squished_Star[:3]+' '+squished_Star[3:]}])
            response_df = response_df.append([{'aliasdis':squished_Star}],ignore_index=True)
        elif first_result != None:
            temp_df = pd.DataFrame([{'aliasdis':first_result}])    
        else: 
            return

        for i in range(response_df.shape[0]):
                response_df.iloc[i] = str(response_df.iloc[i][0]).lower().replace("-","").replace(" ","")       
    
        file_name = temp_df.iloc[0][0]
        temp_df = temp_df.append(response_df,ignore_index=True)
        temp_df = temp_df.T
        for col in temp_df:
            temp_df = temp_df.rename(columns={col:str(col)})
        
        dictionary_df = dictionary_df.append(temp_df,ignore_index=True)
        dictionary_df.to_csv(misc_path +'/Star name dictionary new.csv', index=False)
    else:
        file_name = dictionary_df['0'].iloc[row[0]]
    
    return file_name

#---------------------------------------------------------------------------------------------------------------#
"""DATA VISUALISATION"""
#---------------------------------------------------------------------------------------------------------------#

def graphplotting(toplot):
    global Graphframe
    if toplot[0][0] == 'No data avaiable':
        messagebox.showerror('Error','The Periodogram could not be generated due to an error with the lightcurves.')
        return
    
    f = Figure(figsize=(5,5), dpi=100)
    a = f.add_subplot(111)
    a.set_xlabel(toplot[0][4])
    a.set_ylabel(toplot[0][5])
    a.title.set_text(toplot[0][6])

    for i in range(len(toplot)):
        if toplot[i][2] == 'line':
            a.plot(toplot[i][0],toplot[i][1],linewidth=1,color=toplot[i][3])
        else:
            a.plot(toplot[i][0],toplot[i][1],'x',markersize=10,color=toplot[i][3])
    canvas = FigureCanvasTkAgg(f,master=Graphframe)
    canvas.get_tk_widget().pack(fill='both', expand=True)
    canvas.draw()
    
def graphswitching(w,x,y,z,canplot):
    global Graphframe, SectorExoplanetButtonsframe, OtherControlsframe
    global graphtype, sector, exoplanet, configuration, dfgraph
    if canplot == False:
        return
    if w != graphtype:
        for widget in OtherControlsframe.winfo_children():
            widget.destroy()
        graphtype = w
        configuration = z
        OtherControls(OtherControlsframe)
        
    graphtype, sector, exoplanet, configuration = [w,x,y,z]
    
    for widget in Graphframe.winfo_children():
        widget.destroy()
    for widget in SectorExoplanetButtonsframe.winfo_children():
        widget.destroy()
    
    if len(dfgraph) > 1:
        MakeOptions('Exoplanet ',list(range(1,len(dfgraph)+1)),exoplanet,SectorExoplanetButtonsframe)
    MakeOptions('Sector ',dfgraph[exoplanet]['Sector'].tolist(),sector,SectorExoplanetButtonsframe)
    
    toplot = []
    setup = [widget.winfo_children()[1]['textvariable'] for widget in OtherControlsframe.winfo_children()]
    toplot.append(dfgraph[exoplanet].iloc[sector][graphtype][0])
    
    for i in range(len(setup)):
        if setup[i] == 1:
            toplot.append(dfgraph[exoplanet].iloc[sector][graphtype][i+1])
    
    graphplotting(toplot)
    return 

def MakeTable(data):
    global ResultsTable
    for child in ResultsTable.get_children():
        ResultsTable.delete(child)
        
    ResultsTable["columns"] = list(data.columns)
    ResultsTable["show"] = "headings"
    data_rows = data.to_numpy().tolist()
    for column in ResultsTable["columns"]:
        ResultsTable.heading(column,text=column)
    for row in data_rows:
        ResultsTable.insert("","end",values=row)
    ResultsTable.column('Parameter', width=280)
    ResultsTable.column('Value', width=125)
    ResultsTable.column('Error', width=140)
    
def MakeOptions(text,values,current,frame):    
    Optionsnames = StringVar(frame)
    Options = [text + x for x in [str(i) for i in values]]
    Optionsnames.set(Options[current])
    if text == 'Exoplanet ':
        SelectOptions = OptionMenu(frame, Optionsnames, *Options,command = lambda value: getexoplanet(value,Options))
    elif text == 'Sector ':
        SelectOptions = OptionMenu(frame, Optionsnames, *Options,command = lambda value: getsector(value,Options))
    SelectOptions.pack()
    return

def OtherControls(frame):
    global graphtype, configuration
    peak_tog = {'lc graph':'Show Transits:','lc fold':'Show Drop in Flux','periodogram':'Show Peak'}
    Toggle_Options = [peak_tog[graphtype],'Show Fitting:','Show noisemask:','Show excluded:']
    for i in range(len(configuration)):
        ToggleButton(frame, Toggle_Options[i], configuration[i], graphswitching)
    
def getsector(value,Sectors):
    global dfgraph, graphtype, sector, exoplanet, configuration, canplot
    sector = Sectors.index(value)
    graphswitching(graphtype, sector, exoplanet, configuration, canplot)
    return

def getexoplanet(value,Exoplanets):
    global dfgraph, graphtype, sector, exoplanet, configuration, Results_list, canplot
    default = [['lc graph',[0,0,0,0,0]],['lc fold',[0,0]],['periodogram',[0]]]
    if exoplanet != Exoplanets.index(value):
        configuration = [default[i][1] for i in range(3) if default[i][0]==graphtype][0]
        exoplanet = Exoplanets.index(value)
        MakeTable(Results_list[exoplanet])
    graphswitching(graphtype, 0, exoplanet, configuration, canplot)
    return

def MakeModel():
    global canplot, Final_Results
    if canplot == False:
        return
    yesno = messagebox.askokcancel('3D Model','The model will open in your internet browse,it will also kill the GUI, do you wish to proceed?')
    if yesno == True:
        Make_Model.MainProgram(Final_Results)
    return
        
#---------------------------------------------------------------------------------------------------------------#
"""GUI GLOBAL CODE"""
#---------------------------------------------------------------------------------------------------------------#

def linkpack(topack):
    topack.bind('<Enter>', changecolour)
    topack.bind('<Leave>', revertcolour)
    if str(topack['textvariable']).startswith('https://') == True:
        topack.bind('<Button-1>',lambda e: webbrowser.open_new(str(topack['textvariable'])))
    else:
        topack.bind('<Button-1>',findtutorial)
           
def changecolour(e):
    #This makes the hovered over link appear highlighted#
    e.widget['fg']='blue'
    e.widget['font']="{Arial Nova Cond Light} 14 underline"
    
def revertcolour(e):
    #This sets the colour back to normal once the link isn't hovered over#
    e.widget['fg']='dark blue'
    e.widget['font']="{Arial Nova Cond Light} 14"
    
def findtutorial(e):
    global window, Mainframe
    #This finds the size of the frame and all of its children#
    window.update()
    mycanvas = Mainframe.winfo_children()[0]
    Scroll_frame = mycanvas.winfo_children()[0]
    Frame_size = Scroll_frame.winfo_height()
    scroll_children = Scroll_frame.winfo_children()
    
    #This returns the textvariable of the link which will match with the index of the child#
    if str(e.widget['textvariable']).isnumeric() == True:
        textvar = str(e.widget['textvariable'])
        tut = scroll_children[int(textvar)]
    #This is for linking from another page#
    else:
        PageSwitch(0)
        #This uses the superlink format to find the right text variable#
        textvar = str(e.widget['textvariable']).replace('link','super')
        for i in range(len(scroll_children)):
            #Try/except needed as Frames don't have text variables so throw an error#
            try: 
                if str(scroll_children[i]['textvariable']) == textvar: 
                    tut = scroll_children[i]
                    break
            except:pass
            
    yco = tut.winfo_y()
    #Divides the y co-ordinate of the widget by the height of the frame#
    #And sets the fraction the scrollbar needs to scroll to#
    mycanvas.yview_moveto(yco/Frame_size)

def PageSwitch(Page):
    global currentpage, container_frame
    container_frame.winfo_children()[currentpage].pack_forget()
    currentpage = Page
    
    if Page == 0:
        Mainframe.pack(ipadx=800,ipady=1800,anchor='nw') 
    if Page == 1:
        Page1frame.pack(ipadx=800,ipady=1800,anchor='nw')
    if Page == 2:
        Page2frame.pack(ipadx=800,ipady=1800,anchor='nw')
    if Page == 3:
        Settingsframe.pack(ipadx=800,ipady=1800,anchor='nw')
    return

def txt_to_GUI(frame,Page):
    global childcount, misc_path, Equations
    GUI_text = open(misc_path+'/GUI_text.txt', "r")
    PageType = ['MAIN PAGE',
                'SINGLE TARGET PAGE',
                'MULTIPLE TARGETS PAGE',
                'SETTINGS']
    Frames = [frame]
    paragraph = ''
    childcount = 0
    supercount = 0
    eqcount = 0
    write_mode = False
    for lines in GUI_text.readlines():
        line = str(lines.strip())
        if line == Page: write_mode = True
        elif line in PageType and write_mode == True: break
   
        elif line[:5] == 'mode:' and write_mode == True:
            line = line.replace('mode:','')
            mode = line[:line.index(":")+1]
            input_text = line.replace(mode,'')
            
            #This allows the code to check the mode and change colour, font size, etc#
            if mode == 'set colour:': colour = input_text
            elif mode == 'set size:': font_size = input_text
            elif mode == 'set alignment:': alignment = input_text
            elif mode == 'set frame:': 
                use_frame = Frames[int(input_text)]
                if use_frame == frame: padx = 20
                else: padx = 10
            elif mode == 'set underline:':
                if input_text == 'True': underline = 'underline'
                else: underline = '' 
                
            elif mode == 'create frame:':                    
                #This creates an invisible frame that widgets can be put in#
                new_frame = LabelFrame(use_frame,text=input_text,
                                       font="{} {} {}".format('{Arial Nova Cond Light}',font_size,underline),fg=colour,borderwidth=0)
                new_frame.pack(padx=10,side=alignment,anchor='nw')
                Frames.append(new_frame)
                if use_frame == frame: childcount +=1
                
            elif mode == 'title:':
                title_text = input_text[:input_text.rindex(":")]
                linktype = input_text.replace(title_text+":",'')
                #This creates a title#
                title = Label(use_frame,text=title_text,
                              font="{} {} {}".format('{Arial Nova Cond Light}',font_size,underline),fg=colour)
                packingtime(title, False, padx, 10, alignment,True)
                
                if linktype == 'super': title.config(textvariable = linktype + str(supercount))
                elif linktype == 'not super': title.config(textvariable = linktype + str(childcount))
                if linktype == 'super' or linktype == 'not super':
                    shortcut = Label(Frames[-1],text='- '+title_text,fg=colour,cursor="hand2",textvariable=childcount-1)
                    packingtime(shortcut,True,60,0,'top',False) 
            
            elif mode == 'superlink:':
                shortcut_text = input_text[:input_text.index(":")]
                supervar = input_text.replace(shortcut_text+":",'')
                shortcut = Label(use_frame,text=shortcut_text,fg=colour,cursor="hand2",textvariable=supervar)
                packingtime(shortcut,True,(5,0),0,alignment,False)
                    
        elif write_mode == True:
            #This is just for basic text#
            if line == 'create equation':
                equation = Equations[eqcount]
                eqcount +=1
                childcount += 1
                equation_label = Label(use_frame,image=equation)
                equation_label.pack(anchor='center')
                equation_label.image = equation
            elif line == 'new paragraph':
                paragraph = paragraph + '\n'
            elif line == '':
                vartext = Label(use_frame,text=paragraph,justify='left',wraplengt=1738)
                packingtime(vartext, False, (padx,0), 0, alignment,True)
                paragraph=''
            elif paragraph == '' or paragraph.endswith('\n'): paragraph = paragraph + line 
            else: paragraph = paragraph + ' ' + line
            
def packingtime(topack,link,indent,spacing,side,increase_count):
    global childcount
    #This keeps track of the number of children in the frame so that links do to the right place#
    if increase_count == True: childcount +=1
    #packs the widgets to the specification#
    topack.pack(padx=indent,pady=spacing,anchor='w',side=side)
    if link == True:
        #adds a link binding if it is a link to make it change colour and send you to the correct section#
        linkpack(topack)                
    
#---------------------------------------------------------------------------------------------------------------#
"""GUI BUTTONS CODE"""
#---------------------------------------------------------------------------------------------------------------#  

def ToggleButton(master,title,default,warning):
    global OnOff
    #Custom On/Off button#
    frame = Frame(master)
    #Description of what the button does#
    label = Label(frame,text = title)
    label.pack(side='left')
    #Creates button and sets textvariable either to 0 or 1#
    toggle_btn = Label(frame,image=OnOff[default],textvariable=default)
    toggle_btn.pack(side='left')
    #Pins image to Label#
    toggle_btn.image = OnOff[default]
    #Binds a function to it#
    window.update()
    child = len(master.winfo_children())-1
    toggle_btn.bind('<Button-1>',lambda e: toggle(e,OnOff,default,warning,child))
    frame.pack(padx=20,pady=(5,0),anchor='nw')
    
    return toggle_btn

def toggle(e,OnOff,default,warning,child):
    global graphtype, sector, exoplanet, configuration, canplot
    #Switches between On/Off when pressed#
    #Displays a warning when switch to turn off of default#
    if e.widget['textvariable'] == 1:
        e.widget.config(image=OnOff[0],textvariable=0)
        e.widget.image = OnOff[0]
        if type(warning) != str and warning != None:
            configuration = [0 if i == child else configuration[i] for i in range(len(configuration))]
            graphswitching(graphtype, sector, exoplanet, configuration, canplot)
        elif warning != None and default == 1:
            messagebox.showwarning('Warning',warning)
        
    else:
        e.widget.config(image=OnOff[1],textvariable=1)
        e.widget.image = OnOff[1]
        if type(warning) != str and warning != None:
            configuration = [1 if i == child else configuration[i] for i in range(len(configuration))]
            graphswitching(graphtype, sector, exoplanet, configuration, canplot)
        elif warning != None and default == 0:
            messagebox.showwarning('Warning',warning)
        
def PipelineOptions(frame,gettype):
    global filetype
    #User selects pipeline#
    datatype = StringVar(frame)
    datatype.set(filetype[0])
    #Calls gettype from within desired page class so that it is usable and unique for each page#
    SelectDataType = OptionMenu(frame, datatype, *filetype,command= gettype)
    SelectDataType.pack(side='left',padx=5,ipady=3,pady=5,anchor='nw')
    SelectDataType.config(width=10)
    
    return SelectDataType

def ToggleLock(tolock,lock):
    for widget in tolock:
        if lock == True and type(widget) == Button:
            if widget['text'] == 'Run' or widget['text'] == 'Update' or widget['text'] == 'Cancel Update':
                widget.config(relief='sunken')
            else:
                widget.config(state='disabled')
        elif lock == True:
            widget.config(state='disabled')
        elif lock == False and type(widget) == Button:
            if widget['text'] == 'Run' or widget['text'] == 'Update' or widget['text'] == 'Cancel Update':
                widget.config(relief='raised')
            else:
                widget.config(state='normal')    
        elif lock == False:
            widget.config(state='normal')

#---------------------------------------------------------------------------------------------------------------#
"""GUI PAGE CODE"""
#---------------------------------------------------------------------------------------------------------------#

class MainPage():
    def __init__(self,master):
        #Loads the page on initial startup#
        #Makes the main page scrollable to fit in all the text#
        self.mycanvas = Canvas(master)
        self.mycanvas.pack(side='left',anchor='center',fill='both',expand='yes')
        yscrollbar = ttk.Scrollbar(master, orient='vertical',command=self.mycanvas.yview)
        yscrollbar.pack(side='right', fill='y')
        self.mycanvas.configure(yscrollcommand=yscrollbar.set)
        self.mycanvas.bind('<Configure>', lambda e: self.mycanvas.configure(scrollregion = self.mycanvas.bbox('all')))
        Scroll_frame = Frame(self.mycanvas)
        self.mycanvas.create_window((0,0),window=Scroll_frame,anchor='nw')
        #This binds the Mouse Scrollwheel to the scrollbar#
        self.mycanvas.bind_all("<MouseWheel>", self.Scroll_mousewheel)
        
        #This loads all the text on the Main page from a text file#
        #I put all the text in a text file and read from it to stop it clogging up the code#
        txt_to_GUI(Scroll_frame,'MAIN PAGE')
        
    def Scroll_mousewheel(self, event):
        global currentpage
        #Makes sure the page only scrolls if the user is on the main page#
        if currentpage == 0:
            self.mycanvas.yview_scroll(-1*(int(event.delta/120)), "units")

class SingleTargetPage():
    def __init__(self,master):
        global exoplanet, sector, canplot
        exoplanet = 0
        sector = 0
        self.currenttype = 0 
        txt_to_GUI(master, 'SINGLE TARGET PAGE')
        
        setup_frame = master.winfo_children()[-2]
        self.star_box = Entry(setup_frame, bd=8,font='{Arial Nova Cond Light} 14')
        self.star_box.pack(side='left',padx=(10,0),pady=10,anchor='nw')
        self.Pipelines = PipelineOptions(setup_frame,self.GetType)
        self.run_btn = Button(setup_frame, text='Run',command=lambda: threading.Thread(target=self.SingleTarget).start())
        self.run_btn.pack(side='left',pady=5,anchor='nw')
        self.processing = Label(setup_frame)
        self.processing.pack(side='left',anchor='w',padx=10)
        
        Controls_frame = master.winfo_children()[-1].winfo_children()[0]
        lk_btn = Button(Controls_frame, text ='Lightcurve',command = lambda: graphswitching('lc graph', sector, exoplanet, [0,0,0,0], canplot))
        lk_btn.pack(side='left',anchor='nw')
        flk_btn = Button(Controls_frame, text ='Folded Lightcurve',command = lambda: graphswitching('lc fold', sector, exoplanet, [0,0], canplot))
        flk_btn.pack(side='left',anchor='nw')
        pg_btn = Button(Controls_frame, text ='Periodogram',command = lambda: graphswitching('periodogram', sector, exoplanet, [0], canplot))
        pg_btn.pack(side='left',anchor='nw')
        ThreeD_btn = Button(Controls_frame,text='3D Model',command= lambda: MakeModel())
        ThreeD_btn.pack(side='left',anchor='nw')
        
        Graphcontrols.pack(padx=(20,0),pady=(5,15),side='left',anchor='nw')
        Graphcontrols.pack_propagate(0)
        Graphframe.pack(padx=20,pady=(5,15),side='left',anchor='nw')
        Graphframe.pack_propagate(0)
        Tableframe.pack(padx=(0,20),pady=(5,15),side='left',anchor='nw')
        Tableframe.pack_propagate(0)
        
    def GetType(self,value):
        #Sets the datatype to search for the right type of lightcurve#
        global filetype
        self.currenttype = filetype.index(value)
        return
    
    def SingleTarget(self):
        global dfgraph,runningpage, dictionary_df, firstgraph, canplot, lc_list
        
        if runningpage != 0:
            messagebox.showwarning('Warning','The program is currently running please wait for it to finish.')
            return
        
        ToggleLock([self.run_btn,self.star_box,self.Pipelines],True)
        
        userinput = self.star_box.get()
        file_name = FindName(userinput)
        #Returns if star is not in dictionary#
        if file_name == None: 
            messagebox.showerror('Warning','Could not find star.')
            ToggleLock([self.run_btn,self.star_box,self.Pipelines],False)
            return
                    
        runningpage = 1
        canplot = False
        success = RunProgram(file_name,self.currenttype,self.processing)

        if success == True: firstgraph.set(1)
        runningpage = 0
        self.processing.config(text='')
        ToggleLock([self.run_btn,self.star_box,self.Pipelines],False)
            
            
class MultipleTargetsPage():
    def __init__(self,master):
        self.filename = 0
        self.currenttype = 0
        
        txt_to_GUI(master, 'MULTIPLE TARGETS PAGE')
        choose_file_frame = master.winfo_children()[-5]
        choose_save_loc_frame = master.winfo_children()[-4]
        how_many_frame = master.winfo_children()[-3]
        setup_frame = master.winfo_children()[-2].winfo_children()[0]
        self.loading_frame = master.winfo_children()[-1]
        
        self.choose_file_btn = Button(choose_file_frame, text ='Choose file',command=self.choosefile)
        self.choose_file_btn.pack(padx=12,pady=(20,10),anchor='nw',side='left')
        self.file_frame = LabelFrame(choose_file_frame,height=60,width=800)
        self.file_frame.pack(pady=(20,10),padx=5)
        self.file_frame.pack_propagate(0)
        
        self.choose_save_btn = Button(choose_save_loc_frame, text='Choose save location',command=self.choosesaveloc)
        self.choose_save_btn.pack(padx=12,pady=10,anchor='nw',side='left')
        self.save_frame = LabelFrame(choose_save_loc_frame,height=60,width=700)
        self.save_frame.pack(pady=10,padx=5)
        self.save_frame.pack_propagate(0)
        
        self.how_many = Entry(how_many_frame, bd=8,font='{Arial Nova Cond Light} 14',width=5)
        self.how_many.pack(side='left',padx=(10,0),pady=(6,10),anchor='w')
        
        self.Pipelines = PipelineOptions(setup_frame,self.GetType)
        self.run_btn = Button(setup_frame, text='Run',command = lambda: threading.Thread(target=self.MultipleTargets).start())
        self.run_btn.pack(side='left',pady=5,anchor='nw')
        
        self.progress = ttk.Progressbar(self.loading_frame, orient = 'horizontal', length = 1000, mode = 'determinate')
        self.progress.pack(padx=20,pady=20,anchor='w',side='left')
        self.processing = Label(self.loading_frame)
        self.processing.pack(side='left',anchor='w')
        
        
    def choosefile(self):
        for widget in self.file_frame.winfo_children():
            widget.destroy()
        self.filename = askopenfilename()
        display_text = self.filename.replace(self.filename[:self.filename.rindex("/")]+'/','')
        file = Label(self.file_frame,text=display_text)
        file.pack(padx=5,side='left',anchor='w')
    
    def choosesaveloc(self):
        for widget in self.save_frame.winfo_children():
            widget.destroy()
        self.saveloc = askdirectory()
        save = Label(self.save_frame,text=self.saveloc)
        save.pack(padx=5,side='left',anchor='w')
    
    def GetType(self,value):
        #Sets the datatype to search for the right type of lightcurve#
        global filetype
        self.currenttype = filetype.index(value)
        return
    
    def FormatRepeat(self,repeat,starlist):
        proceed = 'yes'
        try:
            if int(repeat) > len(starlist): 
                repeat = len(starlist)
                proceed = messagebox.askquestion("Max value exceeded","Value exceeds maximum number of stars so all will be analysed. Do you wish to proceed?")
        except ValueError: 
            if repeat.lower() == 'all': repeat = len(starlist)
            else:messagebox.showerror("Warning","Invalid entry.")
        
        return repeat, proceed
    
    def MultipleTargets(self):
        global runningpage, MT_df
        repeat = self.how_many.get()
        
        MT_df = pd.DataFrame(columns = ['Star', 'Star Radius', 'Planet Classification', 'Period (d)','Period Error','Drop in flux','Drop in flux Error','Orbital radius (AU)','Orbital radius Error','Transit time (hrs)','Transit Error',
                                     'Planetary mass (M⊕)','Mass Error1','Mass Error2','Planetary density (g/cm\u00b3)','Density Error1','Density Error2','Planetary radius (R⊕)','Planetary radius Error',
                                     'Planetary velocity (km/s)','Velocity Error','Surface temperature (K)','Surface temperature Error','Surface gravity (m/s\u00b2)', 'Surface gravity Error1','Surface gravity Error2','Within habitable zone'])
        
        if runningpage != 0:
            messagebox.showwarning('Warning','The program is currently running please wait for it to finish.')
            return
        elif self.filename == 0:
            messagebox.showerror("Information","Please select a file.")
        elif len(repeat) == 0:
            messagebox.showerror("Information","Please specify how many stars you want to analyse.")
            
        if self.filename.endswith('.csv'):
            try:
                with open(self.filename, newline='') as f:
                    reader = csv.reader(f)
                    lines = list(reader)
                    
                tic_id_check = [ind for ind in range(len(lines)) if 'tic_id' in lines[ind]]
                ticid_check = [ind for ind in range(len(lines)) if 'ticid' in lines[ind]]
                
                if len(tic_id_check) != 0:
                    dft = pd.read_csv(self.filename, delimiter=',', encoding='latin-1', skiprows=tic_id_check[0])
                    dft_col_list = dft['tic_id'].values.tolist() 
                    starlist = list(dict.fromkeys(dft_col_list))
                    
                elif len(ticid_check) != 0:
                    dft = pd.read_csv(self.filename, delimiter=',', encoding='latin-1', skiprows=ticid_check[0])
                    dft_col_list = dft['ticid'].values.tolist() 
                    dft_col_list = ['TIC '+str(star) for star in dft_col_list]
                    starlist = list(dict.fromkeys(dft_col_list))
                else:
                    messagebox.showerror('Error','Could not find stars, try changing the name of the star column to either tic_id or ticid.')
                    return
            except:
                messagebox.showerror('Error','Could not read file.')
                return
                        
        elif self.filename.endswith('.txt'):
            f = open(self.filename)
            starlist = [x.strip() for x in f.readlines()]
        else:
            messagebox.showerror("Warning","Invalid file type.")
            
        repeat, proceed = self.FormatRepeat(repeat,starlist)
            
        if proceed == 'yes':
            ToggleLock([self.run_btn,self.how_many,self.Pipelines,self.choose_save_btn,self.choose_file_btn], True)
            save_file = datetime.now().strftime(self.saveloc+"/Exoplanet_Detector_Analysis-%Y%m%d%H%M%S.csv")
            runningpage = 2
            
            for i in range(int(repeat)):
                print(i+1,'/',int(repeat))
                self.processing.config(text='Processing...')
                Star = FindName(starlist[i])
                if Star == None: continue
                RunProgram(Star, self.currenttype,self.processing)
                MT_df.to_csv(r""+save_file, index=False)
                percent_done = (i+1)/int(repeat)*100
                self.progress['value'] = percent_done
                window.update_idletasks()
                
            runningpage = 0
            ToggleLock([self.run_btn,self.how_many,self.Pipelines,self.choose_save_btn,self.choose_file_btn], False)
            self.processing.config(text='Done!')

class SettingsPage():
    def __init__(self,master):
        global OnOff
        txt_to_GUI(master,'SETTINGS')
        self.Update_frame = master.winfo_children()[1].winfo_children()[0]
        self.Up_DB_only = ToggleButton(self.Update_frame,'Update existing stars in database:', 1,None)
        self.Up_new = ToggleButton(self.Update_frame,'Search for new stars:', 0,'This can make updating take longer.')
        
        self.Update_btn = Button(self.Update_frame,text='Update',command = lambda: threading.Thread(target=self.Update).start())
        self.Update_btn.pack(side='left',anchor='w',padx=10)
        self.cancel_btn = Button(self.Update_frame,text='Stop Update',command = self.cancelupdate) 
        self.progress = ttk.Progressbar(self.Update_frame, orient = 'horizontal', length = 800, mode = 'determinate')
        self.progress.pack(padx=20,pady=20,anchor='w',side='left')
        self.processing = Label(self.Update_frame)
        self.processing.pack(side='left',anchor='w')
            
    def Update(self):
        global database_path,lc_list,runningpage
        lc_list = []
        conditions = [self.Up_DB_only['textvariable'],self.Up_new['textvariable']]
        if runningpage != 0:
            messagebox.showwarning('Warning','The program is currently running please wait for it to finish.')
            return
        elif conditions == [0,0]:
            messagebox.showwarning('Warning','You have not asked to update anything.')
        else:
            runningpage = 3
            yesno = messagebox.askokcancel('Update database','Are you sure you want to update the Data base, it may take some time.')
            
            if yesno == True and conditions[0] == 1:
                starlist = [file for file in os.listdir(database_path) if os.path.isdir(database_path+'/'+file)==True]
            if yesno == True and conditions[1] == 1:
                new_stars_df = pd.read_csv("https://archive.stsci.edu/missions/tess/catalogs/toi/tois.csv",skiprows=4) 
                new_stars = dict.fromkeys(['TIC '+str(x) for x in new_stars_df['TIC'].tolist()])
                starlist = starlist + [x for x in new_stars if x not in starlist]
            if yesno == True: 
                ToggleLock([self.Update_btn], True)
                self.cancel = 0
                self.cancel_btn.pack(side='left',anchor='w',padx=10)
                for i in range(len(starlist)):
                    if self.cancel == 1:
                        self.cancel_btn.pack_forget()
                        break
                    percent_done = ((i+1)/len(starlist))*100
                    self.progress['value'] = percent_done
                    window.update_idletasks()
                    CheckforNewData(starlist[i], 'ALL',self.processing)
                    
                self.processing.config(text='Done!')
                ToggleLock([self.Update_btn,self.cancel_btn], False)
            runningpage = 0
            
    def cancelupdate(self):
        ToggleLock([self.cancel_btn], True)
        self.processing.config(text='Stoping update...')
        self.cancel = 1
            
#---------------------------------------------------------------------------------------------------------------#
"""STARTUP CODE"""
#---------------------------------------------------------------------------------------------------------------#

class StartUp():
    def __init__(self,master):
        Title = Label(master,text='THE EXOPLANET DETECTOR',font="{Arial Nova Cond Light} 25 underline")
        Title.pack(padx=40,pady=10,anchor='nw')  
        Buttonframe = Frame(master)
        Buttonframe.pack(padx=40,anchor='nw')
        Main_btn = Button(Buttonframe, text ='Home', command= lambda: PageSwitch(0))
        Main_btn.pack(anchor='nw',side='left') 
        Page1_btn = Button(Buttonframe, text ='Single target mode', command= lambda: PageSwitch(1))
        Page1_btn.pack(anchor='nw',side='left')
        Page2_btn = Button(Buttonframe, text ='Multiple target mode', command= lambda: PageSwitch(2))
        Page2_btn.pack(anchor='nw',side='left')
        Settings_btn = Button(Buttonframe, text='Settings',command=lambda: PageSwitch(3))
        Settings_btn.pack(anchor='nw')
        
def firstgraph_callback(self, *args):
    global canplot
    canplot = True
    graphswitching('lc graph', 0, 0, [0,0,0,0], canplot)
        
#All the global variables#
#Works out where the application is saved to find the correct folder# 

# =============================================================================
# try: 
#     os.system("TASKKILL /F /IM ExoplanetDetectorStartup.exe") 
# except:
#     pass
# =============================================================================

application_path = os.path.dirname(os.path.abspath(__file__)) 
misc_path = application_path + '\Misc'
database_path = application_path + '\Exoplanet Light Curves'
filetype = ['FFI',
            'SPOC',
            'TESS-SPOC']
dictionary_df = pd.read_csv(misc_path +'/Star name dictionary.csv')
currentpage = 0
runningpage = 0
graphtype = None
canplot = False
checked_stars =[]

#GUI set up#
window = Tk()
window.state('zoomed')
window.wm_title("The Exoplanet Detector")
defaultFont = font.nametofont("TkDefaultFont")
defaultFont.config(family="Arial Nova Cond Light",size=14)

#Global variable but has to be called after Tk()#
OnOff= [PhotoImage(file=misc_path+'\OFF.png'),PhotoImage(file=misc_path+'\ON.png')]
Equations = [PhotoImage(file=misc_path+'/Equations/'+eq) for eq in os.listdir(str(misc_path+'\Equations'))]
firstgraph = IntVar(window,value=0)
firstgraph.trace("w",firstgraph_callback)
#Calls start up to create general buttons and titles#
StartUp(window)

#Creates the frame which each page will be drawn in, this frame can be used to quickly pack and unpack pages#
container_frame = Frame(window,height=800,width=1800)
container_frame.pack(padx=40,pady=10,anchor='nw')

#Creates the Main page and packs it so that it is the first shown#
Mainframe = LabelFrame(container_frame,text='Home',height=800,width=1800)
Mainframe.pack(ipadx=800,ipady=1800,anchor='nw')

#Creates the other pages but does not pack them#
Page1frame = LabelFrame(container_frame,text='Single target mode',height=800,width=1800)
Graphframe = LabelFrame(Page1frame,height=600,width=870)
Graphcontrols = LabelFrame(Page1frame,height=600,width=330)
SectorExoplanetButtonsframe = Frame(Graphcontrols,height=200)
SectorExoplanetButtonsframe.pack(pady=(20,5))
OtherControlsframe = Frame(Graphcontrols)
OtherControlsframe.pack(anchor='nw')

Tableframe = LabelFrame(Page1frame,height=600,width=550)
ttk.Style().configure('Treeview', rowheight=40)
ttk.Style().configure('Treeview.Heading',font=('Arial Nova Cond Light',14),justify='left',fg='darkblue')
ResultsTable = ttk.Treeview(Tableframe)
ResultsTable.place(relheight=1,relwidth=1)


Page2frame = LabelFrame(container_frame,text='Multiple Targets Mode',height=800,width=1800)

Settingsframe = LabelFrame(container_frame,text='Settings',height=800,width=1800)

#Calls all the page classes to build the widgets for each page#
MainPage(Mainframe)
SingleTargetPage(Page1frame)
MultipleTargetsPage(Page2frame)
SettingsPage(Settingsframe)

if os.path.exists(database_path) == False:
    messagebox.showwarning('Warning','You do not have the database downloaded, without it you will have to download the data for each star. This will cause the program to run much slower. \
A database folder will be created for you in which any stars you analysis will be saved.')
    os.mkdir(database_path)

window.mainloop()
        
        
        
        
        