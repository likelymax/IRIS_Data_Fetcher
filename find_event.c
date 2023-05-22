#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct{
	int day;
	int hour;
	int min;
	int sec;
}ti;

typedef struct{
	int day;
	int hour;
	int min;
	int sec;
	int year;
	int mon;
}tran;

ti findtime(float period){
	ti out;
	out.day = period/3600/24;
	out.hour = (period - out.day * 3600 * 24)/3600;
	out.min = (period - out.day * 3600 * 24 - out.hour * 3600)/60;
	out.sec = period - out.day * 3600 * 24 - out.hour * 3600 - out.min * 60;
	return out;
}

tran transfer(int esec, int emin, int ehour, int eday, int emon, int eyear, ti period){
	int nsec, nmin, nhour, nday, nmon, nyear;
	tran out;
	nsec = esec + period.sec;
	if (nsec > 60) {
		nsec = nsec - 60;
		emin = emin + 1;
	}
	nmin = emin + period.min;
	if ( nmin > 60 ){
		nmin = nmin - 60;
		ehour = ehour + 1;
	}
	nhour = ehour + period.hour;
	if (nhour > 24){
		nhour = nhour - 24;
		eday = eday + 1;
	}
	nday = eday + period.day;
	nmon = emon;
	nyear = eyear;
	if (emon == 2){
		if (eyear == 1996 || eyear == 2000 || eyear == 2004 || eyear == 2008 || eyear == 2012 || eyear == 2016 || eyear == 2020){
			if (nday > 29){
				nday = nday - 29;
				nmon = 3;
			}
		}
	}else if (emon == 1 || emon == 3 || emon == 5 || emon == 7 || emon == 8 || emon == 10){
		if (nday > 31){
			nday = nday - 31;
			nmon = emon + 1;
		}
	}else if (emon == 12){
		if (nday > 31){
			nday = nday - 31;
			nmon = 1;
			nyear = eyear + 1;
		}
	}else if (emon == 4 || emon == 6 || emon == 7 || emon == 9 || emon == 11){
		if (nday > 30){
			nday = nday - 30;
			nmon = emon + 1;
		}
	}
	out.sec = nsec;
	out.min = nmin;
	out.hour = nhour;
	out.day = nday;
	out.mon = nmon;
	out.year = nyear;
	return out;
}

int main(int argc, char **argv){
	int i;
	char event[128], output[128];
	int eyear, emon, eday, ehour, emin, esec;
	int fyear, fmon, fday, fhour, fmin, fsec;
	int length;
	float period, delta;
	ti extra;
	tran endtime;
	FILE *fe, *fo;
	for (i = 0; i < argc; i++){
		if (argv[i][0] == '-'){
			switch (argv[i][1]){
				case 'E':
					strcpy(event, &argv[i][2]);
					fe = fopen(event, "r");
					break;
				case 'T':
					sscanf(&argv[i][2],"%d/%d/%dT%d:%d:%d\n", &eyear, &emon, &eday, &ehour, &emin, &esec);
					break;
				case 'D':
					sscanf(&argv[i][2], "%f/%f\n", &period, &delta);
					break;
				case 'O':
					strcpy(output, &argv[i][2]);
					fo = fopen(output, "w");
					break;			
			}
		}
	}
	int day, year;
	int last;
	int tyear, tmon, tday, thour, tmin, tsec, syear, smon, sday, shour, smin, ssec;
	float mag, evla, evlo, evdp;
	int tim1, tim2, tim3;

	day = esec + emin * 100 + ehour * 10000;
	year= eday + emon * 100 + eyear * 10000; 
	extra = findtime(period);
	endtime = transfer(esec, emin, ehour, eday, emon, eyear, extra);
	fprintf(fo, "%d %d %d %d %d %d\n",endtime.year, endtime.mon, endtime.day, endtime.hour, endtime.min, endtime.sec), 
	fclose(fe);
}