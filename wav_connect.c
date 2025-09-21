#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
// #define M_PI 3.14
typedef struct WAVE_HEADERE
{
	char RIFF[4];                 //RIFF header
	int  total_Len;             //檔案(header+語音)長度 = 36 + dataSize
	char WAVE[4];                 //WAVE header
	
	char FMT[4];                  //fmt header
	int  FMTLen;               
	short int pcm;              
	short int  channels;              //聲道數
	int samplehz;               //取樣點/秒
	int byterate;               //位元速率 = 取樣頻率*位元深度/8
	short int sample_size;            //一個取樣多聲道資料塊大小(bytes)
	short int sample_bits;            //取樣位元深度(m)
	
	
	char DATA[4];                 //data header
	int  DATALen;               //語音資料的大小
}WaveHeader;


void generateWav(FILE *fp,int fs,int m ,double T){


		int sampTimes = fs * T;                    // 次數 = (次數/時間) * 時間
		int dataSize = sampTimes * 1 * (m / 8);    // 大小 = 取樣次數 * 聲道數 * 2byte
		int i = 0;

		short *cosinedata = 0;
		double  P = 0, P0 = 0;                     

		cosinedata = (int16_t*)malloc(sizeof(int16_t)*sampTimes);                        //生成儲存指定cos波的動態陣列
		
		WaveHeader Wh;              //生成空的wav標頭檔以便寫入

		Wh.RIFF[0]='R';                  //RIFF
		Wh.RIFF[1]='I';
		Wh.RIFF[2]='F';
		Wh.RIFF[3]='F';   

		Wh.WAVE[0]='W';                  //WAVE
		Wh.WAVE[1]='A';
		Wh.WAVE[2]='V';
		Wh.WAVE[3]='E';

		Wh.FMT[0]='f';                   //fmt
		Wh.FMT[1]='m';
		Wh.FMT[2]='t';
		Wh.FMT[3]=' ';

		Wh.FMTLen = 16;                  
		Wh.pcm = 1;                   
		Wh.channels = 1;                   //單聲道
		Wh.samplehz = fs;                  //每秒取幾個點
		
		Wh.byterate = fs * 1 * (m / 8);	   //byterate每秒數據字節數= SampleRate * channels * BitsPerSample / 8
		Wh.sample_size = 1 * (m / 8);      //sample_size每個採樣所需的字節數= channels * BitsPerSample / 8 
		          
		Wh.sample_bits = m;                //幾bit             

		Wh.DATA[0]='d';                  //data
		Wh.DATA[1]='a';
		Wh.DATA[2]='t';
		Wh.DATA[3]='a';

		Wh.DATALen = dataSize;             //資料大小
                Wh.total_Len = 36 + dataSize;      //檔案長度 = 36 + dataSize
		fwrite(&Wh,sizeof(Wh),1,fp);                                         //write header
		
		
 
                //生成cos波
                //by formular
		double fi[10]={0,31.25,500,2000,4000,44,220,440,1760,3960};
		double ai[10]={100,2000,1000,500,250,100,2000,1000,500,250};
					
		 for(i = 0; i < sampTimes; i++){
		        int k=0;//for range 0~9
		        int c=0;//check 8000hz or 16000hz
		        if(sampTimes==8000)
		          c=1;//for 8000hz
		        else
		          c=2;//for 16000hz
		          
		        if(i<=c*800)
		          k=0;
		        else if(i<=c*1600)
		          k=1;
		          else if(i<=c*2400)
		          k=2;
		          else if(i<=c*3200)
		          k=3;
		          else if(i<=c*4000)
		          k=4;
		          else if(i<=c*4800)
		          k=5;
		          else if(i<=c*5600)
		          k=6;
		          else if(i<=c*6400)
		          k=7;
		          else if(i<=c*7200)
		          k=8;
		          else
		          k=9;
			double t= (double)i/(double)fs; //t=sampTimes/Hz
			cosinedata[i]=-1*round(ai[k]*cos(2*M_PI*fi[k]*t));                   
			fwrite(&cosinedata[i],sizeof(int16_t),1,fp);
			fflush(fp);
		}
		//free malloc
		free(cosinedata);
		fclose(fp);
		fflush(fp);
}
typedef struct				//定義複數
{
	double real;//real 
	double imag;//imaginary 
}complex;

void rectspect (char wav[],char txt[],int hz,double windowsize,double dftwsz,double frameinterval,int set){

    FILE *file;
    file = fopen(wav,"rb");
	
                WaveHeader wavHead;
		fread(&wavHead,sizeof(wavHead),1,file);
		int tdatasize = wavHead.DATALen;

                //hz=每秒採樣多少點
		int M = hz*frameinterval;                               	 //Number of points in a frame interval
		int N = hz*dftwsz;                                                //Number of points in a DFT analysis
		int P = hz*windowsize;                                                 //Number of points in a analysis window
		int sampTimes = tdatasize/2;                                     //Number of samples
		int fNum = ceil(sampTimes/M);                                //Number of frames
		int countadd = 0;
                int countmul = 0;
		
    	        double **result = (double**)malloc(sizeof(double*)*fNum);   //2-dim array that store result
		for (int z = 0; z < fNum; z++) {
			result[z] = (double*)malloc(sizeof(double)*N);
    	        }
    	        
    	        
    	        int16_t *Data = (int16_t*)malloc(sizeof(int16_t)*sampTimes); //Save data
		fread(Data,tdatasize,1,file);
		
		
	        double *cos_array = (double*)malloc(sizeof(double)*N);
    	        double *sin_array = (double*)malloc(sizeof(double)*N);
		for(int a = 0;a<N;a++){                                             //create cos and sin data
        	cos_array[a] = cos(2*M_PI*(double)a/N);                             // N=DFT point
        	sin_array[a] = sin(2*M_PI*(double)a/N);
    	        }
    	        
                //void spectrogram(char wav[],int fs,double wsz,double fftsz,double frameinterval,char txt[],int set)
                int16_t *rectangular = (int16_t*)malloc(sizeof(int16_t)*sampTimes);       //Save rectangular windowed cosine data
                FILE *txtfile; //Open txt
                txtfile=fopen(txt,"wb");
                
				int a=0;
				while(a<fNum){
				       
				       int mSample=a*M;    
				       complex cacuDFT;
				       cacuDFT.real = 0.0;          //for caculate DFT                                                       
				       cacuDFT.imag = 0.0; 		                                                                
	
			       for(int k=0;k<N;k++){
				                rectangular[k] =0;//intial
				       }                                                       
					for(int q = 0;(q<=N-1)&&(q<=P);q++){

							rectangular[q] = Data[mSample+q]*1;//range 0~N-1 =1
							countmul=countmul+1;//time plus 1
	
					}
                                         
                                        //N=DFT points
					for(int k = 0;k<=N/2;k++){
				
						cacuDFT.real = 0.0;                                                                 
						cacuDFT.imag = 0.0; 
						                                                                //X[m,k] = 20log10|DFT of xm[n]|                                                            
						for(int n = 0;n<N;n++){                                         //用尤拉公式計算dft
							if((n<M)&&(mSample+n<sampTimes)){
								cacuDFT.real = cacuDFT.real+(double)rectangular[n]*cos_array[(k*n)%N]; //complex real
								cacuDFT.imag = cacuDFT.imag-(double)rectangular[n]*cos_array[(k*n)%N]; //complex imagine
								countadd=countadd+2;//add time plus 2
								countmul=countmul+2;//time plus 2
							}
						}
						
						result[a][k]=20*log10(sqrt(pow(cacuDFT.real,2)+pow(cacuDFT.imag,2)));             //Sum up amplititude
						countmul=countmul+1;
						
						fprintf(txtfile,"%lf ",result[a][k]);   //write into txt                                                
						fflush(txtfile);                        //clear cache                                              
					}
					fprintf(txtfile,"\n");
					a++;
				}
				for (int i = 0; i < fNum; i++)        //free 2維陣列 malloc
					free(result[i]);
					
		printf("%s set %d add times= %d mul times= %d\n",wav,set, countadd,countmul);			
	        free(result);                              //free malloc		
		fclose(txtfile);                         //close txt file
                fflush(txtfile);
		free(rectangular);
		free(Data);                                //free malloc
		fclose(file);                                  //close wav file
                fflush(file);     	        
}
void hammingspect (char wav[],char txt[],int hz,double windowsize,double dftwsz,double frameinterval,int set){

    FILE *file;
    file = fopen(wav,"rb");
	
                WaveHeader wavHead;
		fread(&wavHead,sizeof(wavHead),1,file);
		int tdatasize = wavHead.DATALen;

                //hz=每秒採樣多少點
		int M = hz*frameinterval;                               	 //Number of points in a frame interval
		int N = hz*dftwsz;                                                //Number of points in a DFT analysis
		int P = hz*windowsize;                                                 //Number of points in a analysis window
		int sampTimes = tdatasize/2;                                     //Number of samples
		int fNum = ceil(sampTimes/M);                                //Number of frames
		int countadd = 0;
                int countmul = 0;
		
    	        double **result = (double**)malloc(sizeof(double*)*fNum);   //2-dim array that store result
		for (int z = 0; z < fNum; z++) {
			result[z] = (double*)malloc(sizeof(double)*N);
    	        }
    	        
    	        
    	        int16_t *Data = (int16_t*)malloc(sizeof(int16_t)*sampTimes); //Save data
		fread(Data,tdatasize,1,file);
		
		
	        double *cos_array = (double*)malloc(sizeof(double)*N);              //create cos and sin data
    	        double *sin_array = (double*)malloc(sizeof(double)*N);
		for(int a = 0;a<N;a++){                                             
        	cos_array[a] = cos(2*M_PI*(double)a/N);                             // N=DFT point
        	sin_array[a] = sin(2*M_PI*(double)a/N);
    	        }
    	        
                //void spectrogram(char wav[],int fs,double wsz,double fftsz,double frameinterval,char txt[],int set)
                //do hamming
                int16_t *hamming = (int16_t*)malloc(sizeof(int16_t)*sampTimes);       //Save hamming windowed cosine data
                FILE *txtfile; //Open txt
                txtfile=fopen(txt,"wb");
                
				int a=0;
				while(a<fNum){
				       
				       int mSample=a*M;    
				       complex cacuhDFT;
				       cacuhDFT.real = 0.0;   //for caculate DFT                                                                
				       cacuhDFT.imag = 0.0; 		                                                                
				       for(int k=0;k<N;k++){
				                hamming[k] =0;//intial
				       }                                                       
					for(int q = 0;(q<=N-1)&&(q<=P);q++){

							hamming[q] = Data[mSample+q]*(0.54-0.46*cos(2*M_PI*(double)q/(P-1)));//range 0~N-1 =1
							countmul=countmul+1;//time plus 1
							countadd=countadd+1;//time plus 1
	
					}                                         
                                        //N=DFT points
					for(int k = 0;k<=N/2;k++){
				
						cacuhDFT.real = 0.0;                                                                 
						cacuhDFT.imag = 0.0; 
						                                                                //X[m,k] = 20log10|DFT of xm[n]|                                                            
						for(int n = 0;n<N;n++){                                         //用尤拉公式計算dft
							if((n<M)&&(mSample+n<sampTimes)){
								cacuhDFT.real = cacuhDFT.real+(double)hamming[n]*cos_array[(k*n)%N]; //complex real
								cacuhDFT.imag = cacuhDFT.imag-(double)hamming[n]*cos_array[(k*n)%N]; //complex imagine
								countadd=countadd+2;//add time plus 2
								countmul=countmul+2;//time plus 2
							}
						}
						
						result[a][k]=20*log10(sqrt(pow(cacuhDFT.real,2)+pow(cacuhDFT.imag,2)));             //Sum up amplititude
						countmul=countmul+1;
						
						fprintf(txtfile,"%lf ",result[a][k]);   //write into txt                                                
						fflush(txtfile);                        //清除快取                                              
					}
					fprintf(txtfile,"\n");
					a++;
				}
				for (int i = 0; i < fNum; i++)        //free 2維陣列 malloc
					free(result[i]);
					
					
					
		printf("%s set %d add times= %d mul times= %d\n",wav,set, countadd,countmul);			
	        free(result);                              //free malloc		
		fclose(txtfile);                         //close txt file
                fflush(txtfile);
		free(hamming);
		free(Data);                                //free malloc
		fclose(file);                                  //close wav file
                fflush(file);     	        
}

int main(int argc, char *argv[]){

	generateWav(fopen("cos-8k.wav","wb"),8000,16,1);//generate coswav
	generateWav(fopen("cos-16k.wav","wb"),16000,16,1);//generate coswav
	
	

	//total 16 ascii file(4*4)
	rectspect("aeueo-8kHz.wav","aeueo-8kHz{Set1}.txt",8000,0.032,0.032,0.01,1);//set1
	hammingspect("aeueo-8kHz.wav","aeueo-8kHz{Set2}.txt",8000,0.032,0.032,0.01,2);//set2
	rectspect("aeueo-8kHz.wav","aeueo-8kHz{Set3}.txt",8000,0.03,0.032,0.01,3);//set3
	hammingspect("aeueo-8kHz.wav","aeueo-8kHz{Set4}.txt",8000,0.03,0.032,0.01,4);//set4
	
	
	
	rectspect("aeueo-16kHz.wav","aeueo-16kHz{Set1}.txt",16000,0.032,0.032,0.01,1);//set1
	hammingspect("aeueo-16kHz.wav","aeueo-16kHz{Set2}.txt",16000,0.032,0.032,0.01,2);//set2
	rectspect("aeueo-16kHz.wav","aeueo-16kHz{Set3}.txt",16000,0.03,0.032,0.01,3);//set3
	hammingspect("aeueo-16kHz.wav","aeueo-16kHz{Set4}.txt",16000,0.03,0.032,0.01,4);//set4	


	rectspect("cos-8k.wav","cos-8k{set1}.txt",8000,0.032,0.032,0.01,1);//set1
	hammingspect("cos-8k.wav","cos-8k{set2}.txt",8000,0.032,0.032,0.01,2);//set2
	rectspect("cos-8k.wav","cos-8k{set3}.txt",8000,0.03,0.032,0.01,3);//set3
	hammingspect("cos-8k.wav","cos-8k{set4}.txt",8000,0.03,0.032,0.01,4);//set4
	
	rectspect("cos-16k.wav","cos-16k{set1}.txt",16000,0.032,0.032,0.01,1);//set1
	hammingspect("cos-16k.wav","cos-16k{set2}.txt",16000,0.032,0.032,0.01,2);//set2
	rectspect("cos-16k.wav","cos-16k{set3}.txt",16000,0.03,0.032,0.01,3);//set3
	hammingspect("cos-16k.wav","cos-16k{set4}.txt",16000,0.03,0.032,0.01,4);//set4
	

	return 0;
}
