#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

///////////////////////////////////////
// Function: trim - Remove "\n"      //
///////////////////////////////////////

int trim(char *line) {
  int end_pos = strlen(line) - 1;

  if (line[end_pos] == '\n') {
    line[end_pos] = '\0';
    return 1;
  }
  return 0;
}

///////////////////////////////////////////////////////
// Function: get_time_cpu - Get CPU time             //
///////////////////////////////////////////////////////

long get_time_cpu() {
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return ru.ru_utime.tv_sec;
}

///////////////////////////////////////////////////////
// Function: get_time - Get time                     //
///////////////////////////////////////////////////////

long get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec;
}

/////////////////////////////////////////
// Function: count_digit - count digit //
/////////////////////////////////////////

int count_digit(long num) {
  int digit = 1;
  int quotient;

  quotient = int(num / 10);

  while (quotient != 0) {
    digit ++;
    quotient = int(quotient / 10);
  }

  return digit;
}  

//////////////////////////////////////////////////////
// Function: revcomp - convert to reverse complement//
//////////////////////////////////////////////////////

void revcomp(char* str) {
  int i, len;
  char c;

  len = strlen(str);

  for(i=0; i<len/2; i++) {
    c = str[i];
    str[i] = str[len-i-1];
    str[len-i-1] = c;
  }

  for(i=0; i<len; i++) {
    if (str[i] == 'A') {
      str[i] = 'T';
    } else if (str[i] == 'T') {
      str[i] = 'A';
    } else if (str[i] == 'G') {
      str[i] = 'C';
    } else if (str[i] == 'C') {
      str[i] = 'G';
    }
  }
}  