# Time-to-event analysis of fish spawning times

Analysis in the manuscript "Evaluating the potential for prezygotic isolation and hybridization between landlocked and anadromous alewife (Alosa pseudoharengus) following secondary contact." 
Littrell,	K.A.,	Ellis,	D.,	Gephard,	S.,	MacDonald,	A.	D.,	Palkovacs,	E.P.,	Scranton,	K,	and	
Post,	D.M.	(2018)	Evaluating	the	potential	for	pre-zygotic	isolation	and	
hybridization	between	landlocked	and	anadromous	alewife	(Alosa	pseudoharengus)	
following	secondary	contact.	Evolutionary	Applications.	doi-10.1111-eva.12645

Original data and files found at https://datadryad.org/resource/doi:10.5061/dryad.7qt4220

## Data	file

Alewife_Spawning_Data
Description:	Spawning	dates	for	young-of-the-year	alewife	in 5	lakes	in	CT	(two	anadromous	and	three	landlocked). 
Variables measured:
Lake	=	spawning	lake
Form	=	life	history	form,	anadromous	(A)	or	landlocked	(L)
Day	=	Julian	Date
Count	=	number	of	young-of-the-year	alewife	laid	as	eggs	on	that	date (spawning)
Temp	=	lake	temp	(°C),	average	of the	first	three	meters	of	the	epilimnion at	the	
deepest	part	of	the	lake

## Analysis file

Analysis_Spawning.R
Reads	the	local	file	“Alewife_Spawning_Data.csv”	and	maximizes	likelihoods	for	suite	 of	time-to-event	models	with	explanatory	variables:	Form,	Lake,	and	Year	using	 source	code	from	“TimeToEvent_Lik.R”.	Calculates	AIC	values.


