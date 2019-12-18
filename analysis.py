import math
# import scipy.stats as s
import sys

### RUNNING THIS SCRIPT:
## use the following command to run this script:
# python3 analysis.py <log file> 
#
## for example:
# python3 analysis.py myfilteredlog.txt
#
## note that the scipy package is only required to run builtin stats functions
## for error-checking purposes. it is otherwise not needed.

HEADERS = ["timestamp", "_app", "_abtest", "version", "load_time", "click_time", "click_id", "session_id"]
class Results:
    
    def __init__(self, file_path, filter_fun):
        # read data rows from file
        log_file = open(sys.argv[1], "r")
        assert log_file.mode == 'r', "Unable to read file '{}'.".format(file_path)
        raw_results = log_file.readlines()

        # list of dictionaries corresponding to each line of data
        self.log_data = [zip(HEADERS, line[0:-1].split(' ')) for line in raw_results]
        self.log_data = [{x:y for x, y in line} for line in self.log_data]

        # convert all timestamps from string to integer
        for row in self.log_data:
            row["load_time"] = int(row["load_time"])
            row["click_time"] = int(row["click_time"])

        # dict {session id : list of tuples as (load, click)}
        self.clicks_per_session = {}

        # dict {session id : list of load times}
        self.loads_per_session = {}

        # process data for every row that passes the filter function
        for row in self.log_data:
            s_id = "{}:{}".format(row["session_id"], row["version"])

            if (filter_fun(row)):
                if row["click_time"] == 0:
                    # Count as session load
                    prev_loads = self.loads_per_session.get(s_id, [])
                    prev_loads.append(row["load_time"])
                    self.loads_per_session[s_id] = prev_loads
                else:
                    # Count as click event
                    prev_clicks = self.clicks_per_session.get(s_id, [])
                    prev_clicks.append((row["load_time"], row["click_time"]))
                    self.clicks_per_session[s_id] = prev_clicks

        # total num sessions
        self.num_sessions = len(self.loads_per_session)

        # calculate clickthrough rate (proportional)
        self.clickthrough = len(self.clicks_per_session) / self.num_sessions

        # calculate avg time to click (in ms)
        self.time_to_click_vals = self.calc_time_to_click()
        assert (len(self.time_to_click_vals) > 0)
        self.time_to_click = sum(self.time_to_click_vals) / len(self.time_to_click_vals)

        # calculate avg dwell time (in ms)
        self.dwell_time_vals = self.calc_dwell_time()
        assert (len(self.dwell_time_vals) > 0)
        self.dwell_time = sum(self.dwell_time_vals) / len(self.dwell_time_vals)

        # calculate return rate (proportional) (relative to sessions with clicks)
        self.return_total = self.calc_num_return()
        self.return_rate = self.return_total / len(self.clicks_per_session)

        # calculate return rate (proportional) (relative to all sessions)
        self.return_rate_all = self.calc_num_return() / self.num_sessions

    def print_all(self):
        "Prints all metrics"
        print("Total # of Sessions: {}".format(self.num_sessions))
        print("Total # of Returning Sessions: {}".format(self.return_total))
        # print("Total # of Dwell Time Counts: {}".format(len(self.dwell_time_vals))) # sanity check
        print("Clickthrough Rate: {}".format(self.clickthrough))
        print("Average Time to First Click (ms): {}".format(self.time_to_click))
        print("Average Dwell Time (ms): {}".format(self.dwell_time))
        print("Return Rate: {}".format(self.return_rate))

    def calc_time_to_click(self):
        "Returns a list of each session's time to first click"
        time_to_click_vals = []

        # for each session, calculate time between load and first click
        # & append this value to the accumulating list
        for clicks in self.clicks_per_session.values():
            assert (len(clicks) > 0)
            time_to_click_vals.append(clicks[0][1] - clicks[0][0])

        return time_to_click_vals

    def calc_dwell_time(self):
        "Returns a list of each session's dwell time"

        dwell_time_vals = []

        # for each session, if there exists a load after the first click (i.e.,
        # a second load), calculate its dwell time and append to list
        for session, clicks in self.clicks_per_session.items():
            assert (len(clicks) > 0)
            first_click_time = clicks[0][1]
            loads_after = self.get_loads_after(session, first_click_time)
            # if there are any loads after the first click moment, calculate
            # that number as dwell time
            if loads_after != []:
                dwell_time_vals.append(loads_after[0] - first_click_time)

        return dwell_time_vals

    def get_loads_after(self, session_id, time):
        "Returns the number of loads for given session after the given time"
        loads = self.loads_per_session.get(session_id, [])
        return list(filter(lambda x: x > time, loads))

    def calc_num_return(self):
        "Returns the number of unique sessions that left the page and then returned"
        counter = 0

        # for each session, check if a load occurs after any of the clicks (i.e., the first click)
        for session, clicks in self.clicks_per_session.items():
            assert (len(clicks) > 0)
            loads_after = self.get_loads_after(session, clicks[0][1])
            # if loads_after is empty, add one to counter
            counter += (loads_after != []) 

        return counter
                
class Stats():

    def __init__(self, results_all, results_a, results_b):
        self.res_all = results_all
        self.res_a = results_a
        self.res_b = results_b

        # calculate t-test results
        self.time_to_click_res = self.time_to_click_t_test()
        self.dwell_time_res = self.dwell_time_t_test()

        # builtin library t-test for comparison purposes :)
        # print("REAL t-test: {}".format(s.ttest_ind(self.res_a.time_to_click_vals, self.res_b.time_to_click_vals)))
        # print("REAL t-test: {}".format(s.ttest_ind(self.res_a.dwell_time_vals, self.res_b.dwell_time_vals)))

        self.confidence_interval = self.time_to_click_confidence()

    def print_all(self):
        "Print all statistics"
        print("Time To Click T-Test: {}".format(self.time_to_click_res))
        print("Confidence Interval: {}, {}".format(*self.confidence_interval))
        print("Dwell Time T-Test: {}".format(self.dwell_time_res))

    def t_test(self, avg_a, avg_b, vals_a, vals_b):
        "Calculate the t-test statistic for the given values and average"
        diff = self.diff(avg_a, avg_b)
        se = self.standard_err(vals_a, vals_b)
        return (diff / se)

    def diff(self, avg_a, avg_b):
        "Calculate the difference of averages"
        return avg_a - avg_b

    def standard_err(self, vals_a, vals_b):
        "Calculate the standard error between the first and second set of values"
        count_a = len(vals_a)
        count_b = len(vals_b)

        sum_sq_a = sum([x * x for x in vals_a])
        sum_sq_b = sum([x * x for x in vals_b])
        sq_sum_a = sum(vals_a) ** 2
        sq_sum_b = sum(vals_b) ** 2
        
        a = (sum_sq_a - (sq_sum_a / count_a))
        b = (sum_sq_b - (sq_sum_b / count_b))

        x = (a + b) / (count_a + count_b - 2)
        y = (1 / count_a) + (1 / count_b)

        return math.sqrt(x * y)

    def time_to_click_t_test(self):
        "Calculate the t-test statistic for time to first click between A & B"
        return self.t_test(self.res_a.time_to_click,
                            self.res_b.time_to_click, 
                            self.res_a.time_to_click_vals, 
                            self.res_b.time_to_click_vals) 

    def dwell_time_t_test(self):
        "Calculate the t-test statistic for dwell time between A & B"
        return self.t_test(self.res_a.dwell_time,
                            self.res_b.dwell_time, 
                            self.res_a.dwell_time_vals, 
                            self.res_b.dwell_time_vals) 

    def get_confidence_int(self, avg_a, avg_b, vals_a, vals_b, table_val):
        "Calculate the confidence interval for the given values"
        se = self.standard_err(vals_a, vals_b)
        adj_se = table_val * se
        diff = self.diff(avg_a, avg_b)
        
        return (diff - adj_se, diff + adj_se)

    def time_to_click_confidence(self):
        "Calculate confidence interval for time to click"
        TABLE_VAL = 1.69 # provided via table using degrees of freedom
        return self.get_confidence_int(self.res_a.time_to_click,
                                        self.res_b.time_to_click, 
                                        self.res_a.time_to_click_vals, 
                                        self.res_b.time_to_click_vals,
                                        TABLE_VAL) 

        
def main():
    if (len(sys.argv) != 2):
        print("usage: python3 analysis <log file>")
        sys.exit(1)

    # COMPILE RESULTS
    res_all = Results(sys.argv[1], lambda _: True)
    res_a = Results(sys.argv[1], lambda row: row["version"] == "A")
    res_b = Results(sys.argv[1], lambda row: row["version"] == "B")

    stats = Stats(res_all, res_a, res_b)
    
    # PRINT RESULTS
    print("RESULTS FOR ALL")
    res_all.print_all()
    print()

    print("RESULTS FOR A")
    res_a.print_all()
    print()

    print("RESULTS FOR B")
    res_b.print_all()
    print()

    # PRINT STATS
    print("STATS FOR ALL")
    stats.print_all()

    sys.exit(0)

main()