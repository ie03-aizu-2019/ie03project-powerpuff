#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <limits>

using namespace std;

#define INF (1<<21)
#define SIZE 1000
#define TRUE 1
#define FALSE 0

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

struct Cross{
    int a_in,a_out,b_in,b_out;
    //double a_x,a_y,b_x,b_y;
    int a,b;

    double x,y;
    int id;
};

typedef struct{
  double x;
  double y;
  int line1;
  int line2;
  int ID;
}Add_Point;


List list[SIZE];
Result result[SIZE];
Qdata qdata[100];         //0<=Q<=100
Point points[1000];       //2<=N<=1000
Line l[500];              //1<=M<=500
//Store store[100];
Cross cross[100];

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
      //((a-c)^2+(b-d)^2
        double dx2 = ((result[list[i].hen[0]-1].rx - result[list[i].hen[1]-1].rx)
                      *(result[list[i].hen[0]-1].rx - result[list[i].hen[1]-1].rx)
                      +(result[list[i].hen[0]-1].ry - result[list[i].hen[1]-1].ry)
                      *(result[list[i].hen[0]-1].ry - result[list[i].hen[1]-1].ry));

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
    MM[5][7]=0;
    MM[6][1]=0;
    MM[6][8]=0;
    MM[7][5]=0;
    MM[7][9]=0;
    MM[8][6]=0;
    MM[8][2]=0;
    MM[8][9]=sqrt(((result[7].rx - result[8].rx) * (result[7].rx - result[8].rx) +
                   (result[7].ry - result[8].ry) *  (result[7].ry - result[8].ry)));
    MM[9][1]=0;
    MM[9][3]=0;
    MM[9][4]=MM[4][9];
    MM[9][7]=0;
    MM[9][8]=MM[8][9];

}

/*********************************************************/
 // 交差判定(crossing detection)
/*********************************************************/
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

                cross[r].id = r;
                cross[r].a_in = l[i].m1;
                cross[r].a_out = l[i].m2;
                cross[r].b_in = l[j].m1;
                cross[r].b_out = l[j].m2;

                cross[r].a = i;
                cross[r].b = j;

                //  cross[r].a_x = l[i].p1.x;

                result[r+N].rx = cross[r].x = vx;
                result[r+N].ry= cross[r].y = vy;
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

            swap(result[i], result[j]);
        }
        if(result[i].rx == result[j].rx){     //x座標が同じならy座標を参考にする

            swap(result[i], result[j]);
        }
    }
}

/*---------------端点同士からなる線分をダイクストラ用の辺とする--------------*/
for(int i=N+1; i<=N+r; i++){
    for(int j=N+1; j<=N+r; j++){
        list[u].hen[0]=i;
        list[u].hen[1]=j;
        u++;
    }
}
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
        }
        else cout << "NA" << endl;        //始点終点の2データがない
    }
}

  /*------------------------------------------------*/
  /* point:  座標                                   */
  /* vertex: 既存の頂点の数                         */
  /* ID: 座標のID                                   */
  /* return: 受け取った座標idのインデックス         */
  /*------------------------------------------------*/

int Index(Add_Point *point, int vertex, int ID)
{

  for(int i=1; i<=vertex; i++){
    //  cout <<  "(pointID:" << point[i].ID << ", "
    //   << point[i].x << ", "<< point[i].y
    //   << ", ID:" << ID <<") "
    //   << ", vertex: " << vertex << endl;

    if(point[i].ID == ID) return i;
  }
  return -1;
}

  /*------------------------------------------------*/
  /* 最小の⻑さの道をどこに作ればよいか求める       */
  /* point: 座標                                    */
  /* road:  道                                      */
  /* newPoint: 新しい地点                           */
  /* return: 新しい地点から最短の道路網の座標       */
  /*------------------------------------------------*/
Add_Point suggestRoad(Add_Point *point, int road[][2], Add_Point newPoint)
{
  double x1, y1;
  double x2, y2;
  double x3, y3;
  double multi1;
  double multi2;
  int indexP;
  int indexQ;
  double MIN;
  int indexMIN;

  double x[100];
  double y[100];
  double dist[100];

  Add_Point notExist = {-1, -1};
  Add_Point connectPoint;


  /*------------------------------------------------*/
  /* 全ての道に対して内積を求める                   */
  /* 内積が'0'ならある道までの最短距離              */
  /* そのなかで一番距離が短いものを見つける         */
  /*------------------------------------------------*/

  for(int i=1; i<=M; i++){
    indexP = Index(point, N, road[i][0]);
    indexQ = Index(point, N, road[i][1]);

    /* 線分端点'P'の座標 */
    x1 = point[indexP].x;
    y1 = point[indexP].y;
    /* 線分端点'Q'の座標 */
    x2 = point[indexQ].x;
    y2 = point[indexQ].y;

    /* 新しい座標 */
    x3 = newPoint.x;
    y3 = newPoint.y;

    /* 乗算:(x2-x1)^2 */
    multi1 = pow(x2-x1, 2.0);
    multi2 = pow(y2-y1, 2.0);

    x[i] = (1 / (multi1+multi2))*(multi1*x3 + multi2*x1 - (x2-x1)*(y2-y1)*(y1-y3));
    y[i] = (((y2-y1) / (x2-x1)) * (x[i]-x1)) + y1;

    dist[i] = sqrt( fabs((x3-x[i]) * (x3-x[i]) + (y3-y[i]) * (y3-y[i])));
  }
  cout << endl;

  MIN = INF;
  indexMIN = 0;
  for(int i=1; i<=M; i++){
    if((MIN>dist[i]) && (x[i]!=x3) && (y[i])!=y3){
      MIN = dist[i];
      indexMIN = i;
    }
  }
  /*------------------------------------------------*/
  /* 新しい地点からの距離が頂点のほうが近かったら   */
  /* dist[i]を頂点との距離に入れ替える              */
  /*------------------------------------------------*/
  /* double PointToPoint[P];
  for(int i=1; i<N; i++){
    double d = ((x3 - point[i].x)*(x3 - point[i].x))+((y3 - point[i].y)*(y3 - point[i].y));
    PointToPoint[i] = sqrt(d);
    if(PointToPoint[i]<dist[i]  && PointToPoint[i]!=0)
      // PointToPoint[i] = dist[i];
    indexMIN = i;
  }
  */
  connectPoint.x = x[indexMIN];
  connectPoint.y = y[indexMIN];

   cout <<  "(" << connectPoint.x << ", " << connectPoint.y << ") " << endl;

  return connectPoint;
}



//***********************************************************
//main
//***********************************************************
int main(){
  int line[100][2];
  Add_Point point[SIZE];
  Add_Point newPoint[SIZE];
  Add_Point suggestPoint;
    cin >> N >> M >> P >> Q;                   // N,M,P,Qデータ入力

    for(int i=1; i<=N; i++){
        cin >> points[i].x >> points[i].y;     // 座標データ入力
        point[i].x = points[i].x;
        point[i].y = points[i].y;
        point[i].ID = i;
    }

    for(int i = 0; i<N; i++){                  // 座標とその番号をresultに格納
        result[i].rx = points[i+1].x;
        result[i].ry = points[i+1].y;
        result[i].rn = i+1;
    }

    for(int i=1; i<=M; i++){
        cin >> l[i].m1 >> l[i].m2;            // Mデータ入力
        l[i].p1 = points[l[i].m1];
        l[i].p2 = points[l[i].m2];
        line[i][0] = l[i].m1;
        line[i][1] = l[i].m2;
    }
    for(int i=1; i<=P; i++){                   // Pデータ入力
      std::cin >> newPoint[i].x >> newPoint[i].y;
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

    /*---------------------- 最適な道を提案 ------------------------*/
    for(int i=1; i<=P; i++){
      suggestPoint = suggestRoad(point, line, newPoint[i]);
    }
    /*--------------------------------------------------------------*/

    Judgement();     //Qデータ判定


    return 0;
}
