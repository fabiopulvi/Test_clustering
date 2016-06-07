
import java.util.ArrayList;
import java.util.Random;
public class KMeans_Ex4a
{
    private static final int NUM_CLUSTERS = 4;    // Total clusters.
    private static final int TOTAL_DATA = 100;      // Total data points.
    static double imbalance =100;
    static double partition_constraint = (double) (imbalance * (TOTAL_DATA+1) / NUM_CLUSTERS);
    static int max_value=10000000;
    static ArrayList<ArrayList<Double>> dataset = new ArrayList<ArrayList<Double>>();

    private static final double SAMPLES[][] = new double[][] {{1.0, 1.0},
            {1.5, 2.0},
            {3.0, 4.0},
            {5.0, 7.0},
            {3.5, 5.0},
            {4.5, 5.0},
            {3.5, 4.5}};



    private static ArrayList<Data> dataSet = new ArrayList<Data>();
    private static ArrayList<Centroid> centroids = new ArrayList<Centroid>();

    private static void initialize()
    {
        System.out.println("Centroids initialized at:");
        Random r = new Random();
        centroids.add(new Centroid(r.nextInt(max_value) , r.nextInt(max_value))); // lowest set.
        centroids.add(new Centroid(r.nextInt(max_value), r.nextInt(max_value))); // highest set.
        centroids.add(new Centroid(r.nextInt(max_value) , r.nextInt(max_value))); // lowest set.
        centroids.add(new Centroid(r.nextInt(max_value), r.nextInt(max_value))); // highest set.
        centroids.add(new Centroid(1,1)); // lowest set.
        centroids.add(new Centroid(6,6));
        centroids.add(new Centroid(6,1));
        centroids.add(new Centroid(1,6));// highest set.
        System.out.println("     (" + centroids.get(0).X() + ", " + centroids.get(0).Y() + ")");
        System.out.println("     (" + centroids.get(1).X() + ", " + centroids.get(1).Y() + ")");
        System.out.print("\n");


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
                    cluster = i;
                    //System.out.println("point:"+newData.X()+","+newData.Y()+", has been assigned to: "+centroids.get(i).X()+","+centroids.get(i).Y());
                }
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
                System.out.println("the new centroid is in:"+centroids.get(i).X()+","+centroids.get(i).Y()+" centroid has elements: "+centroids.get(i).size);
            }
            centroids.get(i).size(totalInCluster);
        }

        // Now, keep shifting centroids until equilibrium occurs.
        while(isStillMoving)
        {
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
                    System.out.println(+totalX+","+totalInCluster);
                    System.out.println(+totalX+","+totalInCluster);
                    centroids.get(i).X((double) totalX / totalInCluster);
                    centroids.get(i).Y((double) totalY / totalInCluster);
                    centroids.get(i).size(totalInCluster);
                    System.out.println("the new centroid is in:"+centroids.get(i).X()+","+centroids.get(i).Y()+" centroid has elements: "+centroids.get(i).size);

                }
                centroids.get(i).size(totalInCluster);
            }
            boolean equal = true;
            for(int i = 0; i < NUM_CLUSTERS; i++) {
                distance = dist_centroid(centroids_old.get(i), centroids.get(i));
                System.out.println("Now we will compare the old centroids: "+centroids_old.get(i).X()+","+centroids_old.get(i).Y()+" with the new ones: "+centroids.get(i).X()+","+centroids.get(i).Y());

                if (distance!=0) equal=false;
            }
            if (equal==true) {

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
                    System.out.println("     (" + dataSet.get(j).X() + ", " + dataSet.get(j).Y() + ")");
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
        System.out.print("\n");
        return;
    }
}