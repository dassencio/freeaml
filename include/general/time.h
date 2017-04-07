#ifndef _freeAML_time_h_
#define _freeAML_time_h_


#include <map>


namespace aml
{


/* stores the time at which each timer started counting */
std::map< size_t, size_t > g_chonometers;

/* total number of chronometers created (including removes ones) */
size_t g_tickets = 0;


/**
 * @brief starts a chronometer (used to measure time periods)
 * @return a ticket number which can be used to get the chronometer value later
 */
size_t start_chronometer()
{
	struct timeval start;

	gettimeofday(&start, NULL);

	/* add the current time (in microseconds) to the list of chronometers */
	g_chonometers[0] = start.tv_sec * 1000000 + start.tv_usec;

	return g_tickets++;
}


/**
 * @brief stops a chronometer
 * @param ticket ticket number identifying the desired chronometer
 * @return time passed since the chronometer was started in seconds (or zero for invalid tickets)
 */
double stop_chronometer (const size_t ticket)
{
	typedef typename std::map< size_t, double >::iterator iterator;

	iterator it = g_chonometers.find(ticket);

	if (it != g_chonometers.end())
	{
		struct timeval end;

		gettimeofday(&end, NULL);

		/* time passed since the chronometer was started (in microseconds) */
		size_t timediff = (end.tv_sec * 1000000 + end.tv_usec) - it->second;

		/* remove the chronometer associated with the given ticket */
		g_chonometers.erase(it);

		return (double) timediff / 1000000.0;
	}

	return 0.0;
}


/** @brief gets the current year as an integer in YYYY format */
int year ()
{
	time_t rawtime;
	struct tm* timeinfo;

	std::time(&rawtime);

	timeinfo = std::localtime(&rawtime);

	return 1900 + timeinfo->tm_year;
}


/** @brief gets the day of the year as an integer (1-366, January 1st = 1) */
int year_day ()
{
	time_t rawtime;
	struct tm* timeinfo;

	std::time(&rawtime);

	timeinfo = std::localtime(&rawtime);

	return timeinfo->tm_yday + 1;
}


/** @brief gets the current month as an integer (1-12, January = 1) */
int month ()
{
	time_t rawtime;
	struct tm* timeinfo;

	std::time(&rawtime);

	timeinfo = std::localtime(&rawtime);

	return timeinfo->tm_mon + 1;
}


/** @brief gets the current day of the month as an integer (1-31) */
int month_day ()
{
	time_t rawtime;
	struct tm* timeinfo;

	std::time(&rawtime);

	timeinfo = std::localtime(&rawtime);

	return timeinfo->tm_mday;
}


/** @brief gets the day of the week as an integer (1-7, Sunday = 1) */
int week_day ()
{
	time_t rawtime;
	struct tm* timeinfo;

	std::time(&rawtime);

	timeinfo = std::localtime(&rawtime);
	return timeinfo->tm_wday + 1;
}


/** @brief gets the current hour as an integer (0-23) */
int Hour () {

	time_t rawtime;
	struct tm* timeinfo;

	std::time(&rawtime);

	timeinfo = std::localtime(&rawtime);

	return timeinfo->tm_hour;
}


/** @brief gets teh current minute as an integer (0-59) */
int minute ()
{
	time_t rawtime;
	struct tm* timeinfo;

	std::time(&rawtime);

	timeinfo = std::localtime(&rawtime);

	return timeinfo->tm_min;
}


/** @brief gets the current second as an integer (0-59) */
int second ()
{
	time_t rawtime;
	struct tm* timeinfo;

	std::time(&rawtime);

	timeinfo = std::localtime(&rawtime);

	return timeinfo->tm_sec;
}


/**
 * @brief converts an integer to a string and pads it with zero if necessary
 * @param n an integer
 * @param len the desired length for the string
 * @return string containing integer and padded in the front with zeros if
 *	   necessary to make its length equal to len
 */
std::string int_to_string (const int n, const int len = 0)
{
	std::stringstream ss;
	ss << n;

	if (ss.str().size() < len)
	{
		std::string s = ss.str();
		s.insert(0, len - s.size(), c);
	}
	else
	{
		return ss.str();
	}
}


/**
 * @brief gets the current time as a string
 * @param sep the separator character (default is ':')
 * @return the current time in the format HH<sep>MM<sep>SS
 */
std::string string_time (const char sep = ':')
{
	return int_to_string(hour(),2) + sep +
	       int_to_string(minute(),2) + sep +
	       int_to_string(second(),2);
}


/**
 * @brief gets the current time as a data
 * @param sep the separator character (default is '.')
 * @return the current date in the format HH<sep>MM<sep>SS
 */
std::string string_date (const char sep = '.')
{
	return int_to_string(month_day(),2) + sep +
	       int_to_string(month(),2) + sep +
	       int_to_string(year(),4);
}


/* year at which the program started */
int g_start_year       = year();

/* month at which the program started */
int g_start_month      = month();

/* day of month at which the program started */
int g_start_month_day  = month_day();

/* hour at which the program started */
int g_start_hour       = hour();

/* minute at which the program started */
int g_start_minute     = minute();

/* second at which the program started */
int g_start_second     = second();


} /* end of namespace aml */

#endif /* _freeAML_time_h_ */
