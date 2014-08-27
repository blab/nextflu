import datetime

def string_to_date(string):
	(y, m, d) = map(int, string.split('-'))
	return datetime.date(year=y, month=m, day=d)

def year_difference(start_date, end_date):
	start_ord = start_date.toordinal()
	end_ord = end_date.toordinal()
	return (end_ord - start_ord) / 365.25

def numerical_date(date):
	"""Takes a calendar date and a numerical dates in terms of years"""
	start_date = datetime.date(year=date.year, month=1, day=1)
	start_ord = start_date.toordinal()
	end_ord = date.toordinal()
	return date.year + (end_ord - start_ord) / 365.25

def string_to_numerical_date(string):
	date = string_to_date(string)
	return numerical_date(date)

def main():
	date1 = string_to_date('2011-06-12')
	date2 = string_to_date('2013-03-15')
	print year_difference(date1, date2)
	print numerical_date(date1)
	print string_to_numerical_date('2014-02-01')

if __name__ == "__main__":
    main()