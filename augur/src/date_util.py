import datetime

def string_to_date(string):
	(y, m, d) = map(int, string.split('-'))
	return datetime.date(year=y, month=m, day=d)

def year_difference(start_date, end_date):
	start_ord = start_date.toordinal()
	end_ord = end_date.toordinal()
	return (end_ord - start_ord) / 365.25

def year_delta(start_date, years):
	days = round(years * 365.25)
	delta = datetime.timedelta(days=days)
	return start_date + delta

def numerical_date(date, format = '%Y-%m-%d'):
	"""Takes a calendar date and a numerical dates in terms of years"""
	if isinstance(date, basestring):
		try:
			date = datetime.datetime.strptime(date, format).date()
		except:
			date = datetime.datetime.strptime(date, '%Y-%m-%d').date()
	start_date = datetime.date(year=date.year, month=1, day=1)
	start_ord = start_date.toordinal()
	end_ord = date.toordinal()
	return date.year + (end_ord - start_ord) / 365.25

def string_to_numerical_date(string):
	date = string_to_date(string)
	return numerical_date(date)

def date_to_day(date):
	if isinstance(date, datetime.date):
		return date.toordinal()
	elif isinstance(date,basestring):
		return datetime.datetime.strptime(date, '%Y-%m-%d').toordinal()
	elif isinstance(date, int):
		return date
	else:
		print "unknown date format", date
		return np.nan

def main():
	date1 = string_to_date('2011-06-12')
	date2 = string_to_date('2013-03-15')
	print year_difference(date1, date2)
	print numerical_date(date1)
	print string_to_numerical_date('2014-02-01')

if __name__ == "__main__":
	main()