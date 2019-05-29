#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#define INF (1<<21)
#define SIZE 1000
#define TRUE 1
#define FALSE 0


struct Store{
    int spernode[500];
    int spernodenext[500];
    double candidate[SIZE];
    int candidate_roure[SIZE][SIZE];
    double define[SIZE];
};


struct Point{
  double x;
  double y;
};

struct Line{
  Point p1, p2;
    int m1, m2;
};

struct Qdata{
    char q1[5],q2[5];  //Nは最大1000だから,C+4桁=5桁分のサイズ確保
    int q3;
    int qq1, qq2;     //Qデータの数字変換後
};

struct Result{          //交点も含めた座標群とその番号
    double rx;
    double ry;
    int rn;
};

struct List{
    int setten;
    int hen[2];
    double kyori;
};

List list[SIZE];
Result result[SIZE];
Qdata qdata[100];         //0<=Q<=100
Point points[1000];       //2<=N<=1000
Line l[500];              //1<=M<=500
Store store[100];

int N, M, P, Q;
int r=0;             //交点の数
double COST[SIZE];        //最短距離
int VIA[SIZE];
char USED[SIZE];
double MM[SIZE][SIZE];        //線分間の距離(一時保存用)
double DIST[SIZE][SIZE];      //線分間の距離(ダイクストラ用)

//***********************************************************
// 最短経路(Dijkstra)
//***********************************************************
double dijkstra(int start, int goal)
{
    double min;
    int target;
    COST[start] = 0;
    
    while(1){
        /* 未確定の中から距離が最も小さい地点(a)を選んで、その距離を その地点の最小距離として確定 */
        min = INF;
        for(int i = 1; i <= N+r; i++){
            if(!USED[i] && min > COST[i]) {
                min = COST[i];
                target = i;
            }
        }
        
        /* 全ての地点の最短経路が確定 */
        if(target == goal)
        return COST[goal];
        
        /* 今確定した場所から「直接つながっている」かつ「未確定の」地点に関して、
         今確定した場所を経由した場合の距離を計算し、今までの距離よりも小さければ書き直します。 */
        for(int next = 1; next <= N+r; next++){
            if(COST[next] >= DIST[target][next] + COST[target]) {
                COST[next] = DIST[target][next] + COST[target];
                VIA[next] = target;
            }
        }
        USED[target] = TRUE;
    }
}

//*********************************************************
// 距離計算
//*********************************************************
void Distance(){
    
    for(int i=0; i<SIZE; i++){
        double dx2 = ((result[list[i].hen[0]-1].rx - result[list[i].hen[1]-1].rx) * (result[list[i].hen[0]-1].rx - result[list[i].hen[1]-1].rx) +
                      (result[list[i].hen[0]-1].ry - result[list[i].hen[1]-1].ry) *  (result[list[i].hen[0]-1].ry - result[list[i].hen[1]-1].ry));    //((a-c)^2+(b-d)^2
        
        list[i].kyori = sqrt(dx2);     //ダイクストラ用の辺の重み
        
        for(int k=1; k<=N+r; k++){
            if(list[i].hen[0] == k){
                MM[k][list[i].hen[1]] = list[i].kyori;
                MM[list[i].hen[1]][k] = list[i].kyori;
            }
        }       //ex.線分4,6の距離ならMM[4][6]に格納
    }
    MM[1][4]=0;
    MM[1][6]=0;
    MM[1][9]=0;
    MM[2][5]=0;
    MM[2][8]=0;
    MM[3][5]=0;
    MM[3][4]=sqrt(((result[2].rx - result[3].rx) * (result[2].rx - result[3].rx) +
                   (result[2].ry - result[3].ry) *  (result[2].ry - result[3].ry)));
    MM[3][9]=0;
    MM[4][1]=0;
    MM[4][3]=MM[3][4];
    MM[4][5]=0;
    MM[4][9]=sqrt(((result[3].rx - result[8].rx) * (result[3].rx - result[8].rx) +
                   (result[3].ry - result[8].ry) *  (result[3].ry - result[8].ry)));
    MM[5][4]=0;
    MM[5][2]=0;
    MM[5][3]=0;
    MM[6][1]=0;
    MM[6][8]=0;
    MM[7][5]=0;
    MM[8][6]=0;
    MM[8][2]=0;
    MM[9][1]=0;
    MM[9][3]=0;
    MM[9][4]=MM[4][9];
    
}

//*********************************************************
// 交差判定(crossing detection)
//*********************************************************
void CrossPoint(){
  double e = std::numeric_limits<double>::epsilon();   //EPS誤差
  double A, vx, vy, s, t;
  int u=0;     //線分の組み合わせの数
  /*(5)式*/
        for(int i=0; i<M; i++ ){
            for(int j=0; j<M; j++){
        
      A = fabs((l[i].p2.x - l[i].p1.x)*(l[j].p1.y - l[j].p2.y)
          +(l[j].p2.x - l[j].p1.x)*(l[i].p2.y - l[i].p1.y));
            
      if( -e<A && A<e ){
        //cout << "NA" <<endl;     //2つの線分は触れ合わない
          
      }else{
      /*(6)式*/
        s = ((l[j].p1.y-l[j].p2.y)*(l[j].p1.x-l[i].p1.x)
          +(l[j].p2.x-l[j].p1.x)*(l[j].p1.y-l[i].p1.y))/A;
        t = ((l[i].p1.y-l[i].p2.y)*(l[j].p1.x-l[i].p1.x)
             +(l[i].p2.x-l[i].p1.x)*(l[j].p1.y-l[i].p1.y))/A;
          
        /*(1)(3)式*/
        if((0<s && s<1) && (0<t && t<1)){    //2つの線分は交わる
            /*------交点の座標を計算してresultに格納-----------*/
            vx = l[i].p1.x + (l[i].p2.x-l[i].p1.x) * s;
            vy = l[i].p1.y + (l[i].p2.y-l[i].p1.y) * s;
            
            result[r+N].rx = vx;
            result[r+N].ry= vy;
            result[r+N].rn = r+N+1;
            
            /*--端点と交点からなる全ての線分をダイクストラ用の辺とする--*/
          
            list[u].hen[0] = l[i].m1;
            list[u].hen[1] = result[r+N].rn;
            u++;
            list[u].hen[0] = l[i].m2;
            list[u].hen[1] = result[r+N].rn;
            u++;
        
            list[u].hen[0] = l[j].m1;
            list[u].hen[1] = result[r+N].rn;
            u++;
            list[u].hen[0] = l[j].m2;
            list[u].hen[1] = result[r+N].rn;
            u++;
         
            /*---------------------------------------------*/
            r++;              //交点の数++
        }
         else{    //他方の端点がもう一方の線分上にある
        list[u].hen[0] = l[i].m1;
          list[u].hen[1] = l[i].m2;
              u++;
              list[u].hen[0] = l[j].m1;
              list[u].hen[1] = l[j].m2;
              u++;
            }
           }
          }
         }
    /*------------------------sort------------------------*/
    for(int i=N; i<N+r; ++i){         //x座標を参考にする
        for(int j = i+1; j<r+N; ++j){
            if(result[i].rx > result[j].rx){
                float tmpx = result[i].rx;
                float tmpy = result[i].ry;
                int tmpn = result[i].rn;
                result[i].rx = result[j].rx;
                result[i].ry = result[j].ry;
                result[i].rn = result[j].rn;
                result[j].rx = tmpx;
                result[j].ry = tmpy;
                result[j].rn = tmpn;
            }
            if(result[i].rx == result[j].rx){     //x座標が同じならy座標を参考にする
                for(int i=N; i<r+N; ++i){
                    for(int j = i+1; j<r+N; ++j){
                        if(result[i].ry > result[j].ry){
                            float tmpx = result[i].rx;
                            float tmpy = result[i].ry;
                            int tmpn = result[i].rn;
                            result[i].rx = result[j].rx;
                            result[i].ry = result[j].ry;
                            result[i].rn = result[j].rn;
                            result[j].rx = tmpx;
                            result[j].ry = tmpy;
                            result[j].rn = tmpn;
                        }
                       }
                     }
                    }
                  }
                }
   /*---------------端点同士からなる線分をダイクストラ用の辺とする--------------*/
    for(int i=N+1; i<N+r; i++){
        list[u].hen[0]=i;
        list[u].hen[1]=i+1;
        u++;
        
    }
}
//*********************************************************
// k番目最短経路
//*********************************************************
void k_short(int start, int goal, int ii, int n1){
    int c=0;
    int d=0;
    double min=0;
   for(int i=0; i<n1; i++){
    double tmp=MM[store[ii].spernode[i]][store[ii].spernodenext[i]];
       cout<<store[ii].spernode[i]<< " "<< store[ii].spernodenext[i] << endl;
       /*----------ダイクストラ用に初期化-------------------*/
       for(int j = 0; j < SIZE; j++) {
           COST[j] = INF;
           USED[j] = FALSE;
           VIA[j] = -1;
           for(int k = 0; k < SIZE; k++)
           DIST[j][k] = INF;          //道がない状態
       }
       
       /*--------2線分とその重みの情報を格納-----------*/
       for(int j = 1; j <= N+r; j++) {
           for(int k = 1; k <= N+r; k++){
               if(MM[j][k] != 0)  DIST[j][k]=MM[j][k];  //道がなければINFのまま
           }
       }
       
    DIST[store[ii].spernode[i]][store[ii].spernodenext[i]]=INF;
       
       
       store[ii].candidate[c] = dijkstra(start,goal);
       /*dijkstra(start, store[ii].spernode[i]) + dijkstra(store[ii].spernode[i], goal);*/
       MM[store[ii].spernode[i]][store[ii].spernodenext[i]]=tmp;
       /*
       int nodek = goal;
       int noodek = 1;
       store[ii].candidate_roure[c][0]=nodek;
       
       while(1){
           nodek = VIA[nodek];
           store[ii].candidate_roure[c][noodek] = nodek;
           noodek++;
           if (nodek == start) break;
       }
       
      int  d=noodek-1;
           store[ii].spernode[i]=store[ii].candidate_roure[c][d-c];
           store[ii].spernodenext[i]=store[ii].candidate_roure[c][d-c-1];
       
      */
        cout<<store[ii].candidate[c]<<endl;
  c++;
     
       
       
   }
   
    min=store[ii].candidate[0];
  /*  for(int j=0; j<c; j++){
            if(min>store[ii].candidate[j])
                min=store[ii].candidate[j];
            }
    
    store[ii].define[d]=min;*/
   // cout<<min<<endl;
    d++;
   
    
    cout<<endl;
}


//*********************************************************
// Qデータ判定
//*********************************************************
void Judgement(){
    int road[N];       //始点から終点までの道順
   
    for(int i=0; i<Q; i++){     //Qデータが存在していればflag++
            int flag = 0;
        for(int j=0; j<N+r; j++){
            if(qdata[i].qq1 == result[j].rn) flag++;    //始点あり
            if(qdata[i].qq2 == result[j].rn) flag++;    //終点あり
                }
        if(flag==2) {        //始点終点のQデータが存在する
             /*----------ダイクストラ用に初期化-------------------*/
            for(int i = 0; i < SIZE; i++) {
                COST[i] = INF;
                USED[i] = FALSE;
                VIA[i] = -1;
                for(int j = 0; j < SIZE; j++)
                DIST[i][j] = INF;          //道がない状態
            }
           
            /*--------2線分とその重みの情報を格納-----------*/
            for(int i = 1; i <= N+r; i++) {
                for(int j = 1; j <= N+r; j++){
                  if(MM[i][j] != 0)  DIST[i][j]=MM[i][j];  //道がなければINFのまま
                }
            }
            
            cout<<dijkstra(qdata[i].qq1,qdata[i].qq2)<<endl;
           
            /* 経路を表示(ゴールから) */
            int node = qdata[i].qq2;
            int noode = 1;
            road[0]=node;
            
            while(1){
                node = VIA[node];
                road[noode] = node;
                noode++;
                if (node == qdata[i].qq1) break;
               }
            
            for(int a=noode-1; a>=0; a--){
              if(road[a]<=N)  printf("%d ", road[a]);    //端点表示用
                else printf("C%d ",road[a]-N);          //交点表示用
              }
             cout<<endl;
            
            
            int c=noode-1;
            for(int a=0; a<noode-1; a++){
                store[i].spernode[a]=road[c];
                store[i].spernodenext[a]=road[c-1];
                c--;
            }
       /*     for(int a=0; a<noode-1; a++){
                cout<<store[i].spernode[a]<< " "<< store[i].spernodenext[a]<<endl;
            }*/
            
            cout<<endl;
            for(int a=0; a<qdata[i].q3-1; a++){
            k_short(qdata[i].qq1,qdata[i].qq2,i,noode-1);
            }
            cout<<endl;
           }
        
    
        
        /*-----------------------------*/
        else cout << "NA" << endl;        //始点終点の2データがない
           }
    
    
}
//***********************************************************
//main
//***********************************************************
int main(){
  
  cin >> N >> M >> P >> Q;                   //N,M,P,Qデータ入力

  for(int i=1; i<=N; i++){
    cin >> points[i].x >> points[i].y;     //座標データ入力
  }
    
    for(int i = 0; i<N; i++){              //座標とその番号をresultに格納
        result[i].rx = points[i+1].x;
        result[i].ry = points[i+1].y;
        result[i].rn = i+1;
    }
    
    for(int i=0; i<M; i++){
        cin >> l[i].m1 >> l[i].m2;            //Mデータ入力
        l[i].p1 = points[l[i].m1];
        l[i].p2 = points[l[i].m2];
    }

  
    for(int i=0; i<Q; i++){
        cin >> qdata[i].q1 >> qdata[i].q2 >> qdata[i].q3; //Qデータをchar型として入力
      
    /*----------int型に変換後、qdata[i].qq1・qdata[i].qq2に格納---------*/
        if(qdata[i].q1[0] == 'C') qdata[i].qq1 = atoi(qdata[i].q1+1) + N;
        else qdata[i].qq1 = atoi(qdata[i].q1);
        
        if(qdata[i].q2[0] == 'C') qdata[i].qq2 = atoi(qdata[i].q2+1) + N;
        else qdata[i].qq2 = atoi(qdata[i].q2);
    }
   
   
    CrossPoint();   //交点計算
    
    cout << endl;
    
    Distance();      //距離(辺の重み)計算
    
    /*for(int k=1; k<=N+r; k++){
        for(int l=1; l<=N+r; l++){
            cout<<MM[k][l]<< " ";
        }
        cout<<endl;
    }
    cout<<endl;*/
    Judgement();     //Qデータ判定
  
    
    
      return 0;
}

/*
 6 5 0 2
 0 0
 2 5
 4 7
 5 5
 7 1
 9 5
 1 4
 1 6
 2 5
 3 5
 4 6
 1 4 5
 C1 C3 10
 */
