
import java.util.ArrayList;
import java.util.Random;
public class KMeans_Ex4a
{
    private static final int NUM_CLUSTERS = 4;    // Total clusters.
    private static final int TOTAL_DATA = 10000;      // Total data points.
    static double imbalance =5000000;
    static double partition_constraint = (double) (imbalance * (TOTAL_DATA+1) / NUM_CLUSTERS);
    static int max_value=10000000;
    static ArrayList<ArrayList<Double>> dataset = new ArrayList<ArrayList<Double>>();
    static int iteration = 0;
    private static final double SAMPLES[][] = new double[][] {{1.0, 1.0},
            {1.5, 2.0},
            {3.0, 4.0},
            {5.0, 7.0},
            {3.5, 5.0},
            {4.5, 5.0},
            {3.5, 4.5}};
/*
    private int findMedoid(Matrix matrix, Matrix cdata, int cluster, int[] clusterid) {
        double minDistance = Double.MAX_VALUE;
        int medoid = -1;
        for (int row = 0; row < matrix.nRows(); row++) {
            if (clusterid[row] == cluster) {
                double distance = metric.getMetric(matrix, cdata, matrix.getWeights(), row, cluster);
                if (distance < minDistance) {
                    medoid = row;
                    minDistance = distance;
                }
            }
        }
        return medoid;
    }

     * Updates the medoid of one cluster
	 *
	 * @param clusterIndex the index of the cluster
	 * @param clusterMember the indices of the instance members in the cluster
	 * @return true if the medoid of this cluster is changed and update the m_Medoids;
	 * 			otherwise return false
	 */
    /*
private boolean updateMedoid(int clusterIndex, ArrayList<Integer> clusterMember) {
    double BestCost = m_DistanceErrors[clusterIndex];
    int NewMedoid = m_Medoids[clusterIndex];
    int ClusterSize = clusterMember.size();
    for (int i=0; i<ClusterSize; i++) {
        double CurrentCost = 0;
        int CurrentMedoid = clusterMember.get(i);
        for (Integer x : clusterMember) {
            CurrentCost += m_DistanceFunction.distance(
                    m_data.instance(CurrentMedoid), m_data.instance(x));
        }
        if (CurrentCost < BestCost) {
            NewMedoid = CurrentMedoid;
            BestCost = CurrentCost;
        }
    }
    if (NewMedoid == m_Medoids[clusterIndex]) {
        return false;	// Not changed
    }
    else {
        m_Medoids[clusterIndex] = NewMedoid;
        m_DistanceErrors[clusterIndex] = BestCost;
        return true;
    }
}
    */

    private static ArrayList<Data> dataSet = new ArrayList<Data>();
    private static ArrayList<Centroid> centroids = new ArrayList<Centroid>();

    private static void initialize()
    {
        System.out.println("Centroids initialized at:");
        Random r = new Random();
       // centroids.add(new Centroid(r.nextInt(max_value) , r.nextInt(max_value))); // lowest set.
       // centroids.add(new Centroid(r.nextInt(max_value), r.nextInt(max_value))); // highest set.


        for(int i = 0; i < NUM_CLUSTERS; i++)
        {
            centroids.add(new Centroid(r.nextInt(max_value), r.nextInt(max_value)));
        }
        /*
        centroids.add(new Centroid(1,1)); // lowest set.
        centroids.add(new Centroid(6,6));
        centroids.add(new Centroid(6,1));
        centroids.add(new Centroid(1,6));// highest set.
        */
        for(int i = 0; i < NUM_CLUSTERS; i++)
        {
            System.out.println("     (" + centroids.get(i).X() + ", " + centroids.get(i).Y() + ")");
        }



        for (int i=0;i<TOTAL_DATA;i++) {
            ArrayList<Double> point = new ArrayList<Double>();
            double x = (double) r.nextInt(max_value);
            point.add(x);
            double y = (double) r.nextInt(max_value);
            point.add(y);
            dataset.add(point);

           ;
        }
        return;
    }

    private static void kMeanCluster()
    {
        final double bigNumber = Math.pow(10, 10);    // some big number that's sure to be larger than our data range.
        double minimum = bigNumber;
        double maximum=0;
        double distance = 0.0;
        double similarity= 0.0;
        int sampleNumber = 0;
        int cluster = 0;
        boolean isStillMoving = true;
        Data newData = null;
        int size=0;
        int max_iterations=10;

        // Add in new data, one at a time, recalculating centroids with each new one. 
      //{
       //     newData = new Data(dataset.get(sampleNumber).get(0), dataset.get(sampleNumber).get(0));
            while(dataSet.size() < TOTAL_DATA)
            {   newData = new Data(dataset.get(sampleNumber).get(0), dataset.get(sampleNumber).get(0));
               // newData = new Data(SAMPLES[sampleNumber][0], SAMPLES[sampleNumber][1]);

            minimum = bigNumber;
            maximum=0;
            for(int i = 0; i < NUM_CLUSTERS; i++)
            {   size=centroids.get(i).size;
               // System.out.println("this centroid starts has elements: "+centroids.get(i).size);
               // distance = dist(newData, centroids.get(i))* (1 - size / partition_constraint);
                similarity = sim(newData, centroids.get(i))* (1 - size / partition_constraint);
                if(similarity > maximum){
                    maximum = similarity;
                    Random r = new Random();
                    cluster = i;
                    //System.out.println("point:"+newData.X()+","+newData.Y()+", has been assigned to: "+centroids.get(i).X()+","+centroids.get(i).Y());
                }
                Random r = new Random();
                cluster=r.nextInt(NUM_CLUSTERS);
            }
            newData.cluster(cluster);
                dataSet.add(newData);

            sampleNumber++;
        }

        // calculate new centroids.
        for(int i = 0; i < NUM_CLUSTERS; i++)
        {
            int totalX = 0;
            int totalY = 0;
            int totalInCluster = 0;
            for(int j = 0; j < dataSet.size(); j++)
            {
                if(dataSet.get(j).cluster() == i){
                    totalX += dataSet.get(j).X();
                    totalY += dataSet.get(j).Y();
                    totalInCluster++;
                }
            }
            if(totalInCluster > 0){
                centroids.get(i).X(totalX / totalInCluster);
                centroids.get(i).Y(totalY / totalInCluster);
                centroids.get(i).size(totalInCluster);
                //System.out.println("the new centroid is in:"+centroids.get(i).X()+","+centroids.get(i).Y()+" centroid has elements: "+centroids.get(i).size);
            }
            centroids.get(i).size(totalInCluster);
        }
        for(int i = 0; i < NUM_CLUSTERS; i++)
        {
            System.out.println(" centroid:"+i+"    (" + centroids.get(i).X() + ", " + centroids.get(i).Y() + ")");
        }
        // Now, keep shifting centroids until equilibrium occurs.

        while(isStillMoving)
        {   iteration++;
            System.out.println("assigning the points");


            // Assign all data to the new centroids

            System.out.println("Iteration");
            for(int i = 0; i < dataSet.size(); i++)
            {
                Data tempData = dataSet.get(i);
                int old_cluster=tempData.cluster();
                minimum = bigNumber;
                maximum=0;
                for(int j = 0; j < NUM_CLUSTERS; j++)
                {   size = centroids.get(j).size;
                    similarity = sim(tempData, centroids.get(j))* (1 - size / partition_constraint);
                    //System.out.println("raw 112 Compared with: "+centroids.get(j).X()+","+centroids.get(j).Y());
                    if(similarity > maximum){
                        maximum = similarity;
                        cluster = j;
                        //System.out.println("raw 116 point:"+tempData.X()+","+tempData.Y()+", has been assigned to: "+centroids.get(j).X()+","+centroids.get(j).Y());

                    }

                }
                tempData.cluster(cluster);

            }

            //backup of the old centroids
            ArrayList<Centroid> centroids_old = new ArrayList<Centroid>();
            for(int i = 0; i < NUM_CLUSTERS; i++)
            {
                centroids_old.add(new Centroid(centroids.get(i).X(), centroids.get(i).Y()));
            }
            // calculate new centroids.
            for(int i = 0; i < NUM_CLUSTERS; i++)
            {
                double totalX = 0;
                double totalY = 0;
                int totalInCluster = 0;
                for(int j = 0; j < dataSet.size(); j++)
                {
                    if(dataSet.get(j).cluster() == i){
                        totalX += dataSet.get(j).X();
                        totalY += dataSet.get(j).Y();
                        totalInCluster++;
                    }
                }
                if(totalInCluster > 0){
                   // System.out.println(+totalX+","+totalInCluster);
                   // System.out.println(+totalX+","+totalInCluster);
                    centroids.get(i).X((double) totalX / totalInCluster);
                    centroids.get(i).Y((double) totalY / totalInCluster);
                    centroids.get(i).size(totalInCluster);
                    //System.out.println("the new centroid is in:"+centroids.get(i).X()+","+centroids.get(i).Y()+" centroid has elements: "+centroids.get(i).size);
                    //System.out.println("the new centroid is in:"+centroids.get(i).X()+","+centroids.get(i).Y()+" centroid has elements: "+centroids.get(i).size);

                }
                centroids.get(i).size(totalInCluster);
            }
            boolean equal = true;
            for(int i = 0; i < NUM_CLUSTERS; i++)
            {
                System.out.println(" centroid:"+i+"    (" + centroids.get(i).X() + ", " + centroids.get(i).Y() + ") with:"+centroids.get(i).size());
            }
            for(int i = 0; i < NUM_CLUSTERS; i++) {
                distance = dist_centroid(centroids_old.get(i), centroids.get(i));
                //System.out.println("Now we will compare the old centroids: "+centroids_old.get(i).X()+","+centroids_old.get(i).Y()+" with the new ones: "+centroids.get(i).X()+","+centroids.get(i).Y());

                if (distance!=0) equal=false;
            }
            if (equal==true) {

                isStillMoving=false;
            }
            if (iteration==max_iterations) {
                isStillMoving=false;
            }


        }
        return;
    }

    /**
     * // Calculate Euclidean distance.
     * @param d - Data object.
     * @param c - Centroid object.
     * @return - double value.
     */
    private static double dist(Data d, Centroid c)
    {
        return Math.sqrt(Math.pow((c.Y() - d.Y()), 2) + Math.pow((c.X() - d.X()), 2));
    }

    private static double sim(Data d, Centroid c)
    {
        return (double) (1/(1+Math.sqrt(Math.pow((c.Y() - d.Y()), 2) + Math.pow((c.X() - d.X()), 2))));
    }

    private static double dist_centroid(Centroid d, Centroid c)
    {
        return Math.sqrt(Math.pow((c.Y() - d.Y()), 2) + Math.pow((c.X() - d.X()), 2));
    }

    private static class Data
    {
        private double mX = 0;
        private double mY = 0;
        private int mCluster = 0;

        public Data()
        {
            return;
        }

        public Data(double x, double y)
        {
            this.X(x);
            this.Y(y);
            return;
        }

        public void X(double x)
        {
            this.mX = x;
            return;
        }

        public double X()
        {
            return this.mX;
        }

        public void Y(double y)
        {
            this.mY = y;
            return;
        }

        public double Y()
        {
            return this.mY;
        }

        public void cluster(int clusterNumber)
        {
            this.mCluster = clusterNumber;
            return;
        }

        public int cluster()
        {
            return this.mCluster;
        }
    }

    private static class Centroid
    {
        private double mX = 0.0;
        private double mY = 0.0;
        private int size =0;

        public Centroid()
        {
            return;
        }

        public Centroid(double newX, double newY)
        {
            this.mX = newX;
            this.mY = newY;
            return;
        }

        public void X(double newX)
        {
            this.mX = newX;
            return;
        }

        public double X()
        {
            return this.mX;
        }

        public void Y(double newY)
        {
            this.mY = newY;
            return;
        }

        public double Y()
        {
            return this.mY;
        }

        public void size(int newsize)
        {
            this.size = newsize;
            return;
        }

        public double size()
        {
            return this.size;
        }
    }

    public static void main(String[] args)
    {
        initialize();
        kMeanCluster();

        // Print out clustering results.
        for(int i = 0; i < NUM_CLUSTERS; i++)
        {
            System.out.println("Cluster " + i + " includes:");
            for(int j = 0; j < TOTAL_DATA; j++)
            {
                if(dataSet.get(j).cluster() == i){
                    //System.out.println("     (" + dataSet.get(j).X() + ", " + dataSet.get(j).Y() + ")");
                }
            } // j
            System.out.println();
        } // i

        // Print out centroid results.
        System.out.println("Centroids finalized at:");
        for(int i = 0; i < NUM_CLUSTERS; i++)
        {
            System.out.println("     (" + centroids.get(i).X() + ", " + centroids.get(i).Y()+" size:"+centroids.get(i).size );
        }
        System.out.print("iterations: "+iteration+"\n");
        return;
    }
}