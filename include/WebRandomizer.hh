#ifndef _WEBRANDOMIZER_H_
#define _WEBRANDOMIZER_H_
#include<stdio.h>
//by RANDOMIZER
int integer(char a) {
        return((int)a-(int)'0');}
int power(int x,int y)
{
        if(y==0)
        return 1;
        return x*power(x,y-1);
}
int random_int(int arr[],int nos)
{
  //char links[]="wget -nv -q http://67.23.25.127/integers/?num=100^&min=1^&max=100000000^&col=100^&base=10 -O new.txt";
  //char links[]="wget \"http://www.random.org/integers/?num=%HowMuchNumbers%&min=%min%&max=%max%&col=1&base=%base%&format=html&rnd=new\" --no-check-certificate -O new.txt";

  //    system(links);
  int status = system("wget \"http://www.random.org/integers/?num=100&min=1&max=100000000&col=100&base=10&format=html&rnd=new\" --no-check-certificate -O new.txt");

  if(status==0){
    //run as normal
        FILE *files;
        files=fopen("new.txt","r");
        int i=0;char ab;
for(i=0;i<6;i++) fscanf(files,"\n");
for(;;)
     {
                 ab=fgetc(files);
         if(ab=='a'){
                 ab=fgetc(files);
         if(ab=='"'){
                 ab=fgetc(files);
         if(ab=='>'){break;}}}}
        int nos_tp=0,j=0;i=0;int flg=0;
     while(!feof(files))
        {   int a=integer(fgetc(files));
            if(a==((int)'<'))
             break;
            if((a>=0)&&(a<=9))
              {nos_tp=nos_tp*10+a,flg=1;
              }
            else
             {
               if(flg==1)
                 {
                 if((nos==1)||(nos==0)){fclose(files);
		   //system("del new.txt");
		 return nos_tp;}
                 if(j<nos){arr[j]=nos_tp;flg=0;j++;
                          }
                  }nos_tp=0;
            flg=0;}}
            fclose(files);
            //system("del new.txt");
            return (1);
  }
  else{
    cout << "WebRandomizer ERROR: wget returned error code: "
	 << status << endl;
    return(-1);
  }
}
#endif // RANDOMIZER_H_INCLUDED
