import javax.swing.JOptionPane;

public class Jacobi_Seidel {
    public static int n=0;
    static float[][] D;
    static float[][] D1;
    static float[][] R;
    static float[][] matriz;
    static float[] B;
    static float[] RX;
    static float[] X;
    static float[] Xi;
    static float[] bRx;
    static double tol;
    static int C;
    static boolean MD=false;
    public static final int MAX_ITERATIONS = 100;  
    private double[][] M;
   

    public static void main(String[] args) {
        String x=JOptionPane.showInputDialog("Ingrese el numero de cifras significativas: ");
        C=Integer.parseInt(x);
        tol=(0.5)*Math.pow(10, 2-C);
        x=JOptionPane.showInputDialog("Ingrese el tamaño de la matriz: ");
        n=Integer.parseInt(x);
        matriz=new float[n][n+1];
        B=new float[n];
        RX=new float[n];       
        X=new float[n];
        bRx=new float[n];
        Xi=new float[n];

        for(int i=0;i<n;i++){
            x=JOptionPane.showInputDialog("ingresa valores del renglon"+(i+1));
            char[] caracter=x.toCharArray();
            int cont=0;String N="";
            for(int j=0;j<caracter.length;j++){
                if(caracter[j]==45){
                N+="-";
                }else{
                    if((caracter[j]>=48 && caracter[j]<=57) || caracter[j]==46){
                        N+=""+caracter[j]; 
                    }else{
                        matriz[i][cont]=Float.parseFloat(N);
                        cont++;
                        N="";
                    }
                }
                if(j<n){
                X[i]=0;
                }
            }
        }
        x=JOptionPane.showInputDialog("ingresa valores matris b");
            char[] caracter=x.toCharArray();
            int cont=0;String N="";
            for(int j=0;j<caracter.length;j++){
                if(caracter[j]==45){
                N+="-";
                }else{
                    if((caracter[j]>=48 && caracter[j]<=57) || caracter[j]==46){
                        N+=""+caracter[j];                        
                    }else{
                        B[cont]=Float.parseFloat(N);                      
                        matriz[cont][n]=Float.parseFloat(N);
                        cont++;
                        N="";
                    }
                }          
            }
            Mdominante();            
            imprimir3(matriz);
            System.out.println("B:");
            imprimir2(B);    
            Mdominante();
            if(MD==true){
                System.out.println("\nresultado por Jacobi:\n");
                MetodoJacobi();
                System.out.println("\n\nresultado por Gauss-Seidel:\n");
                MetodoGaussSeidel(); 
            }else{
                System.out.println("no es dominante");
            }
            System.out.println("\n\nResultado por Gauss-Jordan:\n");
            Gauss_Jordan();   
    }
    
    static void MetodoJacobi(){
        double suma=0,error=0;
        double tiempo=0,tiempS=0;
        long TInicio, TFin; 
        TInicio = System.currentTimeMillis();
        int cont=0;
        do{          
            jacobi(matriz);        
            RX=multi(R,Xi);         
            bRx=resta(B,RX); 
            Xi=multi(D1,bRx);
            for(int i=0;i<n;i++){
               double yi=Xi[i],y=X[i];
               suma=suma+Math.pow(yi-y, 2);
            }
            error=Math.sqrt(suma);
            suma=0;
            X=Xi;
            cont++;
        }while(error>tol);
        System.out.println(cont+" iteraciones");
        for(int i=0;i<n;i++){
            System.out.println("X["+i+"]= "+Xi[i]);
        }
        //System.out.println("x0= "+CS(Xi[0])+" \nx1= " +CS(Xi[1])+" \nx2= " +CS(Xi[2]));
        TFin = System.currentTimeMillis(); 
        tiempo = TFin - TInicio; 
        tiempS=tiempo/1000;
        tiempo=tiempS/60;
        System.out.println("Tiempo de ejecución en segundos: " + tiempS+"\nTiempo de ejecución en minutos: " + tiempo);     
    }
   
    static void Gauss_Jordan(){
    gauss(matriz);
        float[][] GJ=jordan(matriz);     
        imprimir3( jordan(gauss(GJ)));
    }
    static float[][] gauss(float[][] m){
        float k=0;       
        for(int i=0;i<(n);i++){
            if(m[i][i]!=1 && m[i][i]!=0){
                k=m[i][i];
            for(int K=0;K<=n;K++){                      
                    m[i][K]=m[i][K]/k;
                    }B[i]=B[i]/k;
            }            
            if(m[i][i]==1){
                for(int j=i+1;j<n;j++){
                    k=m[j][i];                   
                    for(int K=0;K<=n;K++){                       
                    m[j][K]=m[j][K]-(m[i][K]*k);
                    }
                }
            }           
        }
        return m;
    }
    
    static float[][] jordan(float[][] m){
        int k=n-1,K=n-1;
        float[][] Mt=new float[n][n+1];
        for(int i=0; i<n;i++){
            for(int j=0;j<n;j++){               
                Mt[i][j]=m[k][K];               
                K--;
            }         
            Mt[i][n]=m[k][n];
            k--;K=n-1;
        }
        return Mt;
    }
    
    static void jacobi(float[][]m){
        D=new float[n][n];
        D1=new float[n][n];
        R=new float[n][n];
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){                   
                if(i==j){                       
                D[i][j]=m[i][j];
                D1[i][j]=1/m[i][j];
                R[i][j]=0;
                }else{
                    R[i][j]=m[i][j];
                    D[i][j]=0;
                    D1[i][j]=0;
                }
            }
        }       
    }         
    
    static float[] resta(float[] b,float[] R){
        float[] res=new float[n];
        for(int i=0;i<n;i++){           
                res[i]=b[i]-R[i];           
        }
        return res;
    }
    
    static float[] multiSeidel(float[][] R,float[] X){
        float[] mul=new float[n];    
        int c=0,C=0;
        float sum=0;
        for(int i=0;i<n;i++){                
            for(int k=0;k<n;k++){                       
                sum=sum+(R[i][k]*X[k]);                   
            }
            mul[i]=sum;                    
            sum=0;                            
        }           
        return mul;
    }
    
    static float[] multi(float[][] R,float[] X){
        float[] mul=new float[n];    
        int c=0,C=0;
        float sum=0;
        for(int i=0;i<n;i++){
            for(int k=0;k<n;k++){
                sum=sum+(R[i][k]*X[k]);
            }
            mul[i]=sum;
            sum=0;
        }
        return mul;
    }
    
    static void imprimir(float[][] m,int M){
        String D="";
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                D+=m[i][j]+" ";
            }
            D+="\n";
        }
        System.out.println(D);
    }
    
    static void imprimir3(float[][] m){
        String D="";
        for(int i=0;i<n;i++){
            for(int j=0;j<n+1;j++){
                D+=m[i][j]+" ";
            }
            D+="\n";
        }
        System.out.println(D);
    }
    
    static void imprimir2(float[] m){
        String D="";
        for(int i=0;i<n;i++){           
                D+=m[i]+" ";           
            D+="\n";
        }
        System.out.println(D);
    }
    
    static void Mdominante(){
        int cont=0;
        double sum=0,a=0;
        for(int i=0;i<n;i++){
            sum=0;
            for(int j=0;j<n;j++){               
                 if(i!=j){                   
                    a=Math.pow(matriz[i][j], 2);
                    a=Math.sqrt(a);                   
                    sum=sum+a;      
                }
            }
            a=Math.pow(matriz[i][i], 2);
            a=Math.sqrt(a);         
            if(a>sum){
                cont++;    
            }
            if(cont==n){
            MD=true;
            }
        }   
    }
        
    public static void MetodoGaussSeidel(){
        double x0,x1,x2,e;
        int cont=0;
        double tiempo=0,tiempS=0;
        long TInicio, TFin; //Variables para determinar el tiempo de ejecución
        TInicio = System.currentTimeMillis();
        float a[][]=new float [n][n+1];
        a=matriz;
        float[] seidel=new float[n];
        for(int i=0;i<n;i++){
        seidel[i]=0;
        }
        x1=0.0; x2=0.0;
        do{ 
            e=seidel[1]; 
            for(int i=0;i<n;i++){
                float sum=0;
                for(int j=0;j<n;j++){
                    if(i!=j){
                        sum=sum-seidel[j]*a[i][j];
                    }
                }
                seidel[i]=(a[i][n]+sum)/a[i][i];
            }
        cont++;
        }while(Math.abs(e-seidel[1])>tol);
        System.out.println(cont+" interacciones");
        for(int i=0;i<n;i++){
            System.out.println("X["+i+"]= "+seidel[i]);
        }
        //System.out.println("x0= "+CS(seidel[0])+" \nx1= " +CS(seidel[1])+" \nx2= " +CS(seidel[2]));
        TFin = System.currentTimeMillis(); 
        tiempo = TFin - TInicio; 
        tiempS=tiempo/1000;
        tiempo=tiempS/60;
        System.out.println("Tiempo de ejecución en segundos: " + tiempS+"\nTiempo de ejecución en minutos: " + tiempo);
        
}
    
    static String CS(double X){
    String D=""+X;
    char[] d=D.toCharArray();
    D="";
    for(int j=0;j<=C;j++){
        D+=d[j];
    }
    return D;
    }
}