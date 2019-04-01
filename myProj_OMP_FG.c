#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include <omp.h>

typedef struct coordinates {
  int num;
  int col;
  int row;
}COORD;

typedef struct sequence {
  int start;
  int end;
  int score;
  int count;
  char *s1;
  char *s2;
}SEQ;

typedef struct parameters{
  char* q;
  char* d;
  int q_cnt;
  int d_cnt;
  FILE* fp;
}PARAM; 

typedef struct performance{
  int total_cells;
  int total_traceback_steps;
  double total_scoreboard_time;
  double total_traceback_time;
}PERF;

void reverse(char *s);
SEQ* calc_2d_table(int n_threads,PERF* performance,int q_len,int d_len,char* q,char * d,int match,int mismatch,int gap);
int calc_score(int match,int mismatch,char c1,char c2);
int calc_cell(int up,int left,int diag,int score,int gap);
SEQ* backtracking(PERF* performance,int **array,int max_q,int max_d,char* q,char* d);
void readPair(PARAM * params,int max_size_q,int min_size_q,int size_d);
double gettime();
void print_performance_stats(PERF* performance,double exec_time,int pairs);

int main(int argc, char *argv[]) {
  int i,j;
  char* name;
  char* path;
  int match;
  int mismatch;
  int n_threads;
  double time1,time2;
  int pairs;
  int min_size_q;
  int max_size_q;
  SEQ* result;
  int size_d;
  char * line;
  int gap;
  for (i=0;i<argc;i++){
    if(strcmp(argv[i],"-name")==0){
      name = malloc(sizeof(char)*strlen(argv[i+1]));
      strcpy(name,argv[i+1]);
    }
    else if(strcmp(argv[i],"-path")==0){
      path = malloc(sizeof(char)*strlen(argv[i+1]));
      strcpy(path,argv[i+1]);
    }
    else if(strcmp(argv[i],"-match")==0){
      match = atoi(argv[i+1]);
    }
    else if(strcmp(argv[i],"-mismatch")==0){
      mismatch = atoi(argv[i+1]);
    }
    else if(strcmp(argv[i],"-gap")==0){
      gap = atoi(argv[i+1]);
    }
    else if(strcmp(argv[i],"-threads")==0){
      n_threads = atoi(argv[i+1]);
    }
  }

  
  printf("\n DEFINED VARIABLES : match = %d, mismatch = %d, gap = %d , name = %s , input path = %s",match,mismatch,gap,name,path);
  FILE* fp = fopen ("Report_ID.txt","w");;
  time1 = gettime();
  PERF* performance = (PERF*)malloc(sizeof(PERF));
  performance->total_cells=0;
  performance->total_traceback_time=0;
  performance->total_scoreboard_time=0;
  performance->total_traceback_steps=0;
  PARAM * params = (PARAM*)malloc(sizeof(PARAM));
  params->fp = fopen(path, "r"); 
  params->q = NULL;
  params->d = NULL;
  fscanf(params->fp,"%*s %d\n",&pairs);
  fscanf(params->fp,"%*s %d\n",&min_size_q);
  fscanf(params->fp,"%*s %d\n",&max_size_q);
  fscanf(params->fp,"%*s %d\n",&size_d);
  fflush(stdin);
  printf("\nStarting procedure with params --> pairs : %d || max_size_q : %d || min_size_q : %d || size_d : %d",pairs,max_size_q,min_size_q,size_d);
  for(i=0;i<pairs;i++){
    if(params->q !=NULL){ free(params->q); }
    if(params->q !=NULL){ free(params->d); }
    params->q = NULL;
    params->d = NULL;
    params->q_cnt = 0;
    params->d_cnt = 0;
    readPair(params,max_size_q,min_size_q,size_d);
    result = calc_2d_table(n_threads,performance,params->q_cnt,params->d_cnt,params->q,params->d,match,mismatch,gap);
    for(j=0;j<result[0].count;j++){
      fprintf(fp,"\nQ:\t%s",params->q);
      fprintf(fp,"\nD:\t%s",params->d);
      fprintf(fp,"\n MATCH %d [Score:%d, Start:%d, End:%d]",j+1,result[j].score,result[j].start,result[j].end);
      fprintf(fp,"\n\tQ:%s",result[j].s1);
      fprintf(fp,"\n\tD:%s",result[j].s2);
      fprintf(fp,"\n");
      free(result[j].s1);
      free(result[j].s2);
    }
    free(result);     
  }
  fclose(fp);
  fclose(params->fp);
  time2 = gettime();
  print_performance_stats(performance,time2-time1,pairs);
  return 0;
}

void print_performance_stats(PERF* performance,double exec_time,int pairs){
  printf("\n======================================");
    printf("\nPROGRAM FINISHED ");
    printf("\nA) PAIRS : %d",pairs);
    printf("\nB) CALCULATED CELLS : %d",performance->total_cells);
    printf("\nC) TOTAL TRACEBACK STEPS : %d",performance->total_traceback_steps);
    printf("\nD) TOTAL EXECUTION TIME : %f",exec_time);
    printf("\nE) TOTAL SCOREBOARD FILLING TIME : %f",performance->total_scoreboard_time);
    printf("\nF) TOTAL TRACEBACK TIME : %f",performance->total_traceback_time);
    printf("\nG) CUPS based on execution_time : %f",performance->total_cells/exec_time);
    printf("\nH) CUPS based on scoreboard_time: %f",performance->total_cells/performance->total_scoreboard_time);
    printf("\n======================================\n");
}

void readPair(PARAM * params,int max_size_q,int min_size_q,int size_d){
  int QorD=0;
  char c;
  params->q = (char*)malloc(sizeof(char)*(max_size_q+1));
  params->d = (char*)malloc(sizeof(char)*(size_d+1));
  while ((c = fgetc(params->fp)) != EOF)
    {
      if(c=='Q'){
        if(QorD==1){break;}
        QorD =0;
    } 
    else if(c=='D'){
      QorD=1;
    }
    else{
      if (isalpha(c)||isdigit(c)) {
        if(QorD==0){
          params->q[params->q_cnt]=c;
          (params->q_cnt)++;
        }
        else{
          params->d[params->d_cnt]=c;
          (params->d_cnt)++;
        }     
      }
    }
  }
  params->q[params->q_cnt]='\0';
  params->d[params->d_cnt]='\0';
  char* tmp;
  tmp = (char*)realloc(params->q,sizeof(char)*(params->q_cnt+1));
  params->q = tmp;
}

SEQ* calc_2d_table(int n_threads,PERF* performance,int q_len,int d_len,char* q,char * d,int match,int mismatch,int gap){
  int **array = malloc(sizeof(int*) * (q_len+1));
  int score = 0;
  int currentRow = 0;
  int currentColumn = 0;
  double time1,time2;
  int *curr_ptr;
  int i,j;

  for(i=0;i<= q_len;i++){
    array[i]=malloc(sizeof(int)*(d_len+1));
  }
  for(i=0;i<=d_len;i++){
    array[0][i]=0;
  }
  for(i=0;i<=q_len;i++){
    array[i][0]=0;
  }
  time1 = gettime();

  printf("\nFilling elements of the scoreboard");
  omp_set_num_threads(n_threads);
  #pragma omp parallel default(none) shared(array,q_len,d_len,match,mismatch,gap,d,q) private(i,j,currentRow,currentColumn,score) 
  {
    
    for(i = 0; i <= (q_len+d_len+1); i++) {
      if(q_len>i){
        j=i;
        currentColumn=0;
      }
      else{
        j=q_len;
        currentColumn=i-q_len;
      }
      #pragma omp for
      for(currentRow = j; currentRow >= 0; currentRow--) {
        if(currentColumn > d_len){continue;}
        if(currentColumn>=1 && currentRow>=1){
          score = calc_score(match,mismatch,q[currentRow-1],d[currentColumn-1]);
          array[currentRow][currentColumn] = calc_cell(array[currentRow][currentColumn-1],array[currentRow-1][currentColumn],array[currentRow-1][currentColumn-1],score,gap);
          
          currentColumn++;
        }
      }
    }
  }
  performance->total_cells = (q_len-1)*(d_len-1);
  time2 = gettime();
  performance->total_scoreboard_time += time2-time1;
  time1 = gettime();
  SEQ* sequences = backtracking(performance,array,q_len,d_len,q,d);
  time2 = gettime();
  performance->total_traceback_time += time2-time1;
  
  for (i = 0; i <= q_len; i++)
  {
    curr_ptr = array[i];
    free(curr_ptr);
  }
  return sequences;
}

SEQ* backtracking (PERF* performance,int **array,int max_q,int max_d,char* q,char* d){
  printf("\nstarting backtracking!");
  int max =0;
  COORD * max_array =NULL;
  int max_array_num = 0;
  int tmp = -1;
  int tc;
  int tr;
  int i,j;
  SEQ * sequence_array;
  char* s1 = NULL;
  char* s1_t = NULL;
  char* s2 = NULL;
  char* s2_t = NULL;
  char t1[2],t2[2];
  int cnt=0;
  int start,end;                                      
  for(i=0;i<=max_q;i++){
    for (j=0;j<=max_d;j++){
      if(array[i][j]>max){
        max = array[i][j];
        max_array_num = 0;
        if(max_array != NULL){free(max_array);max_array= NULL;}
        COORD * tmp_array; 
        tmp_array = realloc(max_array,sizeof(COORD)*(++max_array_num));
        max_array = tmp_array;
        max_array[max_array_num-1].col = j;
        max_array[max_array_num-1].row = i; 
        max_array[max_array_num-1].num = array[i][j];
        //printf("\nrow:%d col:%d max_num:%d num:%d",i,j,max_array_num,max_array[max_array_num]);
             
      }
      else if(array[i][j]==max && max!=0){
        COORD * tmp_array = realloc(max_array,sizeof(COORD)*(++max_array_num));                                                                   
        max_array = tmp_array;
        max_array[max_array_num-1].col = j;
        max_array[max_array_num-1].row = i; 
        max_array[max_array_num-1].num = array[i][j];
        //printf("\nrow:%d col:%d max_num:%d",i,j,max_array_num);
      }
    }
  }
  sequence_array =  (SEQ*)malloc(sizeof(SEQ)*(max_array_num));
  //#pragma omp parallel for default(none) shared(max_array,array,sequence_array,q,d,performance,max_array_num) private(tc,tr,t1,t2,tmp,cnt,start,end,s1,s2,s1_t,s2_t,i)
  for(i=0;i<max_array_num;i++){
    tc = max_array[i].col;
    tr = max_array[i].row;
    end = tc;  
    tmp=-1;
    cnt =0;
    while(tmp!=0){
      performance->total_traceback_steps++;
      if(array[tr-1][tc-1]>=array[tr][tc-1] && array[tr-1][tc-1]>=array[tr-1][tc]){
        tmp= array[tr-1][tc-1];
        t1[0] = q[tr-1];
        t1[1] = '\0';
        t2[0] = d[tc-1];
        t2[1] = '\0';
        tc= tc-1;
        tr =tr-1;
        cnt++;
      }
      else if(array[tr][tc-1]>=array[tr-1][tc-1] && array[tr][tc-1]>=array[tr-1][tc]){
        tmp = array[tr][tc-1];
          t1[0] = '-';
          t1[1] = '\0';
          t2[0] =  d[tc-1];
          t2[1] = '\0';
          tc= tc-1;
          tr =tr;
          cnt++;
      }
      else if(array[tr-1][tc]>=array[tr-1][tc-1] && array[tr-1][tc]>=array[tr][tc-1]){
        tmp = array[tr-1][tc];
          t1[0] = q[tr-1];
          t1[1] = '\0';
          t2[0] = '-'; 
          t2[1] = '\0';
          tc= tc;
          tr =tr-1;
          cnt++;
      }
      s1_t = (char*)realloc(s1,sizeof(char)*(strlen(t1)+(cnt)));
      s1 = s1_t;
      s2_t = (char*)realloc(s2,sizeof(char)*(strlen(t2)+(cnt)));
      s2 = s2_t;
      if(cnt==1){
        strcpy(s1,t1);
        strcpy(s2,t2);
      }
      else{
        strcat(s1,t1);
        strcat(s2,t2);
      }
      start = tc;
    }
    sequence_array[i].end = end;
    sequence_array[i].start = start;
    sequence_array[i].score = max_array[0].num;
    sequence_array[i].count = max_array_num;
    sequence_array[i].s1=(char*)malloc((strlen(s1)+1)*sizeof(char));
    strcpy(sequence_array[i].s1,s1);
    reverse(sequence_array[i].s1);
    sequence_array[i].s2=(char*)malloc((strlen(s2)+1)*sizeof(char));
    strcpy(sequence_array[i].s2,s2);
    reverse(sequence_array[i].s2);
    free(s1);
    free(s2);
    s1=NULL;s2=NULL;s1_t=NULL;s2_t=NULL;
  }
  free(max_array);
  return sequence_array;
}

double gettime()
{
struct timeval ttime;
gettimeofday(&ttime, NULL);
return ttime.tv_sec+ttime.tv_usec * 0.000001;
}

void reverse(char *s)
{
   int length, c;
   char *begin, *end, temp;
   length = strlen(s);
   begin  = s;
   end    = s;
   for (c = 0; c < length - 1; c++)
      end++;
   for (c = 0; c < length/2; c++)
   {        
      temp   = *end;
      *end   = *begin;
      *begin = temp;
      begin++;
      end--;
   }
}

int calc_score(int match,int mismatch,char c1,char c2){
  if(c1==c2){
      return match;
  }
  else{
    return mismatch;
  }
}

int calc_cell(int up,int left,int diag,int score,int gap){
  int n1,n2,n3;
  int final_num = 0;
  n1 = up + gap;
  n2 = left + gap;
  n3 = score +diag;
  if( n1>=n2 && n1>=n3 ){final_num = n1;}
  if( n2>=n1 && n2>=n3 ){final_num = n2;}
  if( n3>=n1 && n3>=n2 ){final_num = n3;}
  if(final_num <0){return 0;}else{return final_num;}
}
