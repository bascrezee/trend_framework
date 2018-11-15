import cdstoolbox as ct

@ct.application(title='Extract a time series and download as NetCDF')
@ct.input.dropdown('variable', values=[
    '2m_temperature','volumetric_soil_water_layer_1'
])

@ct.input.dropdown('gridpoint_region',values=['Antarctic','Central Europe'])
@ct.output.dataarray()

def retrieve_time_series(
    variable,
    gridpoint_region
):
    """
    Application main steps:
    
    - retrieve a variable over a defined time range
    - select a location, defined by a region name
    - compute the daily average
    - show the result as a timeseries on an interactive chart
    
    """
    gridpoint_dictionary = {
        'Antarctic' : (59,-75),
        'Central Europe' : (11,50),
    }
    
    longitude,latitude = gridpoint_dictionary[gridpoint_region]
    
    # Time range
    data = ct.catalogue.retrieve(
        'reanalysis-era5-single-levels',
        {
            'variable': variable,
            'grid': ['3', '3'],
            'product_type': 'reanalysis',
            'year': 
                list(range(2000,2019)),
            'month': [
                '01', '02', '03', '04', '05', '06',
                '07', '08', '09', '10', '11', '12'
            ],
            'day': [
                '01', '02', '03', '04', '05', '06',
                '07', '08', '09', '10', '11', '12',
                '13', '14', '15', '16', '17', '18',
                '19', '20', '21', '22', '23', '24',
                '25', '26', '27', '28', '29', '30',
                '31'
            ],
            'time': ['00:00', '06:00', '12:00', '18:00'],
        }
    )    
    
    # Location selection
    
    # Extract the closest point to selected lon/lat (no interpolation). 
    # If wrong number is set for latitude, the closest available one is chosen:
    # e.g. if lat = 4000 -> lat = 90.
    # If wrong number is set for longitude, first a wrap in [-180, 180] is made,
    # then the closest one present is chosen:
    # e.g. if lon = 200 -> lat = -160.
    
    
    data_sel = ct.geo.extract_point(data, lon=longitude, lat=latitude)
    
    # Daily mean on selection
    data_daily = ct.climate.daily_mean(data_sel)
    
    
    return data_daily

