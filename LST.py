def JD(date, UT):
    '''Compute the Julian day give date and time at the Universal Time
        
        Parameters
        ---------------
        year : Int
        month: Int
        day :  Int
        UT (Universal time): float in hours
        
        returns
        ------------
        Day: float
        
        ------------
        #ref: http://cds.cern.ch/record/1634847/files 978-3-319-01442-5_BookBackMatter.pdf
        The formula for JD is valid for the years 1901_2099
        '''
    year, month, day  = date
    if month <= 2:
        year = year - 1
        month = month + 12
    
    Jd=367 * year - (int(7 * (year+int((month+9.)/12.)))/4) + int((275 * month)/9.0) + day + 1721013.5 + UT/24.
    return Jd


def gsmt(JD):
    """Compute the mean sidereal time at Greenwhich given the Julian day
        
        Parameters
        -------------
        JD : -float
        -Julian day
        
        Returns
        --------------
        Time: in hours
        """
    julian_century = 36525 #days
    t = (JD - 2451545.0)/julian_century
    #Greenhich sidereal time at Jan 1 of the year of interest (theta0)
    theta0 = 280.46061837 + (360.98564736629 *(JD  - 2451545.0)) +( 0.00387933 * t * t)-(t*t*t)/38710000.0 #degrees
    
    return ((theta0 % 360.0) * 12.0)/180.0



def LST(JD,longitude):
    """Compute the Local Sidereal Time
        
        when East - add longitude
        when West - substract longitude
        
        Parameters
        -------------------------------------
        longitude: [hrs]
        
        returns:
        LST in hour fraction
        """
    long = rad2hr(longitude)
    gsmst = gsmt(JD)
    lst_HRS=gsmst - long
    
    return lst_HRS


def getHA(ra, lst, units = 'hours'):
    """
        compute the hour angle given the lst [rad] and Right Ascension [rad]
        
        Parameters
        ---------------
        ra : right ascession [radians]
        lst : local sidereal time [radians]
        
        returns:
        --------------
        ha : hour angle [radians]
        """
    
    try: ra, lst
    except NameError:
        raise NameError('The right acession and local sidereal time must be given.')
    
    if isinstance(ra, (tuple,list)):
        ra = np.array(ra)
    
    if units == 'hours':
        lst = hr2rad(lst)

    ha = lst - ra
    if ha > 2 * np.pi: ha -= 2 * np.pi

    return ha


