import math
import re
from datetime import date

# --------------------- Copyright Block ----------------------
# praytimes.py: Prayer Times Calculator (ver 2.3)
# Copyright (C) 2007-2011 PrayTimes.org
# Python Code: Saleem Shafi, Hamid Zarrabi-Zadeh
# Original js Code: Hamid Zarrabi-Zadeh
# License: GNU LGPL v3.0
# TERMS OF USE:
# 	Permission is granted to use this code, with or
# 	without modification, in any website or application
# 	provided that credit is given to the original work
# 	with a link back to PrayTimes.org.
# This program is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY.
# PLEASE DO NOT REMOVE THIS COPYRIGHT BLOCK.

# ----------------------- PrayTimes Class ------------------------


class PrayTimes:
    # ------------------------ Constants --------------------------
    TIME_NAMES = {
        'imsak': 'Imsak', 'fajr': 'Fajr', 'sunrise': 'Sunrise',
        'dhuhr': 'Dhuhr', 'asr': 'Asr', 'sunset': 'Sunset',
        'maghrib': 'Maghrib', 'isha': 'Isha', 'midnight': 'Midnight'
    }

    METHODS = {
        'MWL': {'name': 'Muslim World League', 'params': {'fajr': 18, 'isha': 17}},
        'ISNA': {'name': 'Islamic Society of North America (ISNA)', 'params': {'fajr': 15, 'isha': 15}},
        'Egypt': {'name': 'Egyptian General Authority of Survey', 'params': {'fajr': 19.5, 'isha': 17.5}},
        'Makkah': {'name': 'Umm Al-Qura University, Makkah', 'params': {'fajr': 18.5, 'isha': '90 min'}},
        'Karachi': {'name': 'University of Islamic Sciences, Karachi', 'params': {'fajr': 18, 'isha': 18}},
        'Tehran': {'name': 'Institute of Geophysics, University of Tehran',
                   'params': {'fajr': 17.7, 'isha': 14, 'maghrib': 4.5, 'midnight': 'Jafari'}},
        'Jafari': {'name': 'Shia Ithna-Ashari, Leva Institute, Qum',
                   'params': {'fajr': 16, 'isha': 14, 'maghrib': 4, 'midnight': 'Jafari'}}
    }

    DEFAULT_PARAMS = {'maghrib': '0 min', 'midnight': 'Standard'}

    # ---------------------- Default Settings --------------------
    DEFAULT_SETTINGS = {
        "imsak": '10 min', "dhuhr": '0 min', "asr": 'Standard', "highLats": 'NightMiddle'
    }

    TIME_SUFFIXES = [' am', ' pm']
    INVALID_TIME = '-----'
    NUM_ITERATIONS = 1

    # ---------------------- Initialization -----------------------
    def __init__(self, method="MWL"):
        self.settings = self.DEFAULT_SETTINGS.copy()
        self.offset = {name: 0 for name in self.TIME_NAMES}
        self.calc_method = method if method in self.METHODS else 'MWL'

        for m_name, config in self.METHODS.items():
            for param_name, value in self.DEFAULT_PARAMS.items():
                if param_name not in config['params'] or config['params'][param_name] is None:
                    config['params'][param_name] = value

        self.settings.update(self.METHODS[self.calc_method]['params'])

        self.lat = 0
        self.lng = 0
        self.elv = 0
        self.time_zone = 0
        self.jDate = 0
        self.time_format = '24h'

    # -------------------- Interface Functions --------------------
    def set_method(self, method):
        if method in self.METHODS:
            self.adjust(self.METHODS[method]['params'])
            self.calc_method = method

    def adjust(self, params):
        self.settings.update(params)

    def tune(self, time_offsets):
        self.offset.update(time_offsets)

    def get_method(self):
        return self.calc_method

    def get_settings(self):
        return self.settings

    def get_offsets(self):
        return self.offset

    def get_defaults(self):
        return self.METHODS

    def get_times(self, date_input, coords, timezone, dst=0, time_format=None):
        self.lat, self.lng = coords[0], coords[1]
        self.elv = coords[2] if len(coords) > 2 else 0
        if time_format is not None:
            self.time_format = time_format

        year, month, day = date_input if isinstance(date_input, tuple) else (date_input.year, date_input.month, date_input.day)
        self.time_zone = timezone + (1 if dst else 0)
        self.jDate = self._julian(year, month, day) - self.lng / (15 * 24.0)
        return self._compute_times()

    def get_formatted_time(self, time, time_format, suffixes=None):
        if math.isnan(time):
            return self.INVALID_TIME
        if time_format == 'Float':
            return time

        suffixes = suffixes if suffixes is not None else self.TIME_SUFFIXES
        time = self._fix_hour(time + 0.5 / 60)  # add 0.5 minutes to round
        hours = math.floor(time)
        minutes = math.floor((time - hours) * 60)
        suffix = suffixes[0 if hours < 12 else 1] if time_format == '12h' else ''

        if time_format == "24h":
            formatted_time = f"{int(hours):02d}:{int(minutes):02d}"
        else:
            formatted_time = f"{int((hours + 11) % 12 + 1)}:{int(minutes):02d}"
        return formatted_time + suffix

    # ---------------------- Calculation Functions -----------------------
    def _mid_day(self, time):
        eqt = self._sun_position(self.jDate + time)[1]
        return self._fix_hour(12 - eqt)

    def _sun_angle_time(self, angle, time, direction=None):
        try:
            decl = self._sun_position(self.jDate + time)[0]
            noon = self._mid_day(time)
            t = 1 / 15.0 * self._arccos((-self._sin(angle) - self._sin(decl) * self._sin(self.lat)) /
                                        (self._cos(decl) * self._cos(self.lat)))
            return noon + (-t if direction == 'ccw' else t)
        except ValueError:
            return float('nan')

    def _asr_time(self, factor, time):
        decl = self._sun_position(self.jDate + time)[0]
        angle = -self._arccot(factor + self._tan(abs(self.lat - decl)))
        return self._sun_angle_time(angle, time)

    def _sun_position(self, jd):
        D = jd - 2451545.0
        g = self._fix_angle(357.529 + 0.98560028 * D)
        q = self._fix_angle(280.459 + 0.98564736 * D)
        L = self._fix_angle(q + 1.915 * self._sin(g) + 0.020 * self._sin(2 * g))

        e = 23.439 - 0.00000036 * D

        RA = self._arctan2(self._cos(e) * self._sin(L), self._cos(L)) / 15.0
        eqt = q / 15.0 - self._fix_hour(RA)
        decl = self._arcsin(self._sin(e) * self._sin(L))

        return decl, eqt

    def _julian(self, year, month, day):
        if month <= 2:
            year -= 1
            month += 12
        A = math.floor(year / 100)
        B = 2 - A + math.floor(A / 4)
        return math.floor(365.25 * (year + 4716)) + math.floor(30.6001 * (month + 1)) + day + B - 1524.5

    # ---------------------- Compute Prayer Times -----------------------
    def _compute_prayer_times(self, times):
        times = self._day_portion(times)
        params = self.settings

        imsak = self._sun_angle_time(self._eval(params['imsak']), times['imsak'], 'ccw')
        fajr = self._sun_angle_time(self._eval(params['fajr']), times['fajr'], 'ccw')
        sunrise = self._sun_angle_time(self._rise_set_angle(self.elv), times['sunrise'], 'ccw')
        dhuhr = self._mid_day(times['dhuhr'])
        asr = self._asr_time(self._asr_factor(params['asr']), times['asr'])
        sunset = self._sun_angle_time(self._rise_set_angle(self.elv), times['sunset'])
        maghrib = self._sun_angle_time(self._eval(params['maghrib']), times['maghrib'])
        isha = self._sun_angle_time(self._eval(params['isha']), times['isha'])

        return {
            'imsak': imsak, 'fajr': fajr, 'sunrise': sunrise, 'dhuhr': dhuhr,
            'asr': asr, 'sunset': sunset, 'maghrib': maghrib, 'isha': isha
        }

    def _compute_times(self):
        times = {
            'imsak': 5, 'fajr': 5, 'sunrise': 6, 'dhuhr': 12,
            'asr': 13, 'sunset': 18, 'maghrib': 18, 'isha': 18
        }

        for _ in range(self.NUM_ITERATIONS):
            times = self._compute_prayer_times(times)

        times = self._adjust_times(times)

        if self.settings['midnight'] == 'Jafari':
            times['midnight'] = times['sunset'] + self._time_diff(times['sunset'], times['fajr']) / 2
        else:
            times['midnight'] = times['sunset'] + self._time_diff(times['sunset'], times['sunrise']) / 2

        times = self._tune_times(times)
        return self._modify_formats(times)

    def _adjust_times(self, times):
        params = self.settings
        tz_adjust = self.time_zone - self.lng / 15.0

        for t in times:
            times[t] += tz_adjust

        if params['highLats'] != 'None':
            times = self._adjust_high_lats(times)

        if self._is_min(params['imsak']):
            times['imsak'] = times['fajr'] - self._eval(params['imsak']) / 60.0

        if self._is_min(params['maghrib']):
            times['maghrib'] = times['sunset'] - self._eval(params['maghrib']) / 60.0

        if self._is_min(params['isha']):
            times['isha'] = times['maghrib'] - self._eval(params['isha']) / 60.0

        times['dhuhr'] += self._eval(params['dhuhr']) / 60.0

        return times

    def _asr_factor(self, asr_param):
        methods = {'Standard': 1, 'Hanafi': 2}
        return methods.get(asr_param, self._eval(asr_param))

    def _rise_set_angle(self, elevation=0):
        elevation = 0 if elevation is None else elevation
        return 0.833 + 0.0347 * math.sqrt(elevation)  # an approximation

    def _tune_times(self, times):
        for name in times:
            times[name] += self.offset[name] / 60.0
        return times

    def _modify_formats(self, times):
        for name in times:
            times[name] = self.get_formatted_time(times[name], self.time_format)
        return times

    def _adjust_high_lats(self, times):
        params = self.settings
        night_time = self._time_diff(times['sunset'], times['sunrise'])

        times['imsak'] = self._adjust_hl_time(
            times['imsak'], times['sunrise'], self._eval(params['imsak']), night_time, 'ccw')
        times['fajr'] = self._adjust_hl_time(
            times['fajr'], times['sunrise'], self._eval(params['fajr']), night_time, 'ccw')
        times['isha'] = self._adjust_hl_time(
            times['isha'], times['sunset'], self._eval(params['isha']), night_time)
        times['maghrib'] = self._adjust_hl_time(
            times['maghrib'], times['sunset'], self._eval(params['maghrib']), night_time)
        return times

    def _adjust_hl_time(self, time, base, angle, night, direction=None):
        portion = self._night_portion(angle, night)
        diff = self._time_diff(time, base) if direction == 'ccw' else self._time_diff(base, time)
        if math.isnan(time) or diff > portion:
            time = base + (-portion if direction == 'ccw' else portion)
        return time

    def _night_portion(self, angle, night):
        method = self.settings['highLats']
        portion = 1 / 2.0  # midnight
        if method == 'AngleBased':
            portion = 1 / 60.0 * angle
        elif method == 'OneSeventh':
            portion = 1 / 7.0
        return portion * night

    def _day_portion(self, times):
        return {k: v / 24.0 for k, v in times.items()}

    # ---------------------- Misc Functions -----------------------
    def _time_diff(self, time1, time2):
        return self._fix_hour(time2 - time1)

    def _eval(self, st):
        val = re.split(r'[^0-9.+-]', str(st), 1)[0]
        return float(val) if val else 0.0

    def _is_min(self, arg):
        return isinstance(arg, str) and 'min' in arg

    # ----------------- Degree-Based Math Functions -------------------
    @staticmethod
    def _sin(d): return math.sin(math.radians(d))

    @staticmethod
    def _cos(d): return math.cos(math.radians(d))

    @staticmethod
    def _tan(d): return math.tan(math.radians(d))

    @staticmethod
    def _arcsin(x): return math.degrees(math.asin(x))

    @staticmethod
    def _arccos(x): return math.degrees(math.acos(x))

    @staticmethod
    def _arctan(x): return math.degrees(math.atan(x))

    @staticmethod
    def _arccot(x): return math.degrees(math.atan(1.0 / x))

    @staticmethod
    def _arctan2(y, x): return math.degrees(math.atan2(y, x))

    def _fix_angle(self, angle): return self._fix(angle, 360.0)

    def _fix_hour(self, hour): return self._fix(hour, 24.0)

    @staticmethod
    def _fix(a, mode):
        if math.isnan(a):
            return a
        a -= mode * (math.floor(a / mode))
        return a + mode if a < 0 else a


# ---------------------- prayTimes Object -----------------------
pray_times_calculator = PrayTimes()

# -------------------------- Parameter --------------------------
CALC_METHOD = 'Karachi'
PLACE = 'Dhaka'
LATITUDE = 23.7104
LONGITUDE = 90.40744
GMT = +6
DST = 0  # Daylight Saving Time
TIME_FORMAT = '12h'

# sample code to run in standalone mode only
if __name__ == "__main__":
    print(f'Prayer Times for today in {PLACE}\n{date.today()}\n')
    times = pray_times_calculator.get_times(
        date.today(), (LATITUDE, LONGITUDE), GMT, DST, TIME_FORMAT)

    for i in ['Fajr', 'Sunrise', 'Dhuhr', 'Asr', 'Maghrib', 'Isha', 'Midnight']:
        print(f"{i}: {times[i.lower()]}")
