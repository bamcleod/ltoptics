import CESAPI.connection
import CESAPI.command
import CESAPI.packet
import CESAPI.refract
import CESAPI.test


class LaserTracker():
    def __init__(self, host,logger):
        self.host = host
        self.logger = logger
        #self.connect()
        
    def connect(self):
        self.reconnect(self)
        
    def reconnect(self):
        if hasattr(self, 'connection'):
            print('Disconnecting from laser tracker at ',self.host)
            self.connection.disconnect()
        else:
            print('Initial connection')
            self.connection = CESAPI.connection.Connection()
            self.ri_algorithm = CESAPI.refract.AlgorithmFactory().refractionIndexAlgorithm(CESAPI.refract.RI_ALG_Leica)
        print('Connecting to laser tracker at ',self.host)
        self.stream = self.connection.connect(host=self.host)
        self.command = CESAPI.command.CommandSync(self.connection)
        #print(self.stream, self.command)
        
    def initialize(self):
        global running
        running = False
        manualiof = False
        forceinit = True
        command = self.command
        units = CESAPI.packet.SystemUnitsDataT()
        units.lenUnitType = CESAPI.packet.ES_LU_Millimeter  # ES_LengthUnit
        # units.angUnitType = ES_AU_Radian  # ES_AngleUnit
        # units.tempUnitType = ES_TU_Celsius  # ES_TemperatureUnit
        # units.pressUnitType = ES_PU_Mbar  # ES_PressureUnit
        # units.humUnitType = ES_HU_RH  # ES_HumidityUnit
        self.logger.debug('Setting units...')
        command.SetUnits(units)

        status = command.GetSystemStatus()
        self.logger.debug('Tracker Processor Status: {}'.format(status.trackerProcessorStatus))
        if forceinit or status.trackerProcessorStatus != CESAPI.packet.ES_TPS_Initialized:  # ES_TrackerProcessorStatus
            self.logger.debug('Initializing...')
            command.Initialize()

        self.logger.debug('setting measurement mode...')
        command.SetMeasurementMode(CESAPI.packet.ES_MM_Stationary)  # ES_MeasMode (only choice for AT4xx)

        self.logger.debug('setting stationary mode parameters...')
        mode_params = CESAPI.packet.StationaryModeDataT()
        mode_params.lMeasTime = 1000  # 1 second
        command.SetStationaryModeParams(mode_params)

        self.logger.debug('setting coordinate system type to Right-Handed Rectangular...')
        command.SetCoordinateSystemType(CESAPI.packet.ES_CS_RHR)  # one of ES_CoordinateSystemType enumerated in CESAPI/packet.py
    #    self.logger.debug('setting coordinate system type to spherical SCC... ')
    #    command.SetCoordinateSystemType(CESAPI.packet.ES_CS_SCC)  # spherical coords

        self.logger.debug('setting system settings...')
        settings = CESAPI.packet.SystemSettingsDataT()
        # one of ES_WeatherMonitorStatus
        if manualiof:
            settings.weatherMonitorStatus = CESAPI.packet.ES_WMS_ReadOnly
        else:
            settings.weatherMonitorStatus = CESAPI.packet.ES_WMS_ReadAndCalculateRefractions
        settings.bApplyStationOrientationParams = int(1)
        settings.bKeepLastPosition = int(1)
        settings.bSendUnsolicitedMessages = int(1)
        settings.bSendReflectorPositionData = int(0)
        settings.bTryMeasurementMode = int(0)
        settings.bHasNivel = int(1)
        settings.bHasVideoCamera = int(1)
        command.SetSystemSettings(settings)

        #reset()

    def measure(self):
        CESAPI.refract.SetRefractionIndex(self.command, self.ri_algorithm)
        self.measurement =  self.command.StartMeasurement()
        return ([self.measurement.dVal1, self.measurement.dVal2, self.measurement.dVal3])

    def goto_position(self, position):
        self.command.GoPosition(True, position[0], position[1], position[2])
    
