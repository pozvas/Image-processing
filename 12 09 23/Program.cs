// See https://aka.ms/new-console-template for more information

using System.ComponentModel;
using System.Drawing;
using System.Drawing.Imaging;
using static System.Net.Mime.MediaTypeNames;

namespace image
{
    internal class Program
    {
        static void Main(string[] args)
        {
            /*Image aboba = new Image("C:\\Users\\Василий\\Pictures\\3.jpg", "C:\\Users\\Василий\\Pictures\\3.1.jpg");
            Console.WriteLine(aboba.MSE(0, aboba.image2.Width, 0, aboba.image1.Height));
            Console.WriteLine(aboba.MSE_Part(10, 10));
            Console.WriteLine(aboba.YIK(0, aboba.image2.Width, 0, aboba.image1.Height));
            Console.WriteLine(aboba.YIK_Part(10, 10));

            Console.WriteLine("\n");

            Image aboba3 = new Image("C:\\Users\\Василий\\Pictures\\5.png", "C:\\Users\\Василий\\Pictures\\5.png");
            Console.WriteLine(aboba3.MSE(0 , aboba3.image2.Width, 0, aboba3.image1.Height));
            Console.WriteLine(aboba3.MSE_Part(10, 10));
            Console.WriteLine(aboba3.YIK(0, aboba3.image2.Width, 0, aboba3.image1.Height));
            Console.WriteLine(aboba3.YIK_Part(10, 10));

            Console.WriteLine("\n");
            
            CompareImage aboba2 = new CompareImage("C:\\Users\\Василий\\Pictures\\5.1.png", "C:\\Users\\Василий\\Pictures\\5.1.png");
            Console.WriteLine(aboba2.MSE(0, aboba2.image2.Width, 0, aboba2.image1.Height));
            Console.WriteLine(aboba2.MSE_Part(10, 10));
            Console.WriteLine(aboba2.YIK(0, aboba2.image2.Width, 0, aboba2.image1.Height));
            Console.WriteLine(aboba2.YIK_Part(10, 10));
            */

            Noise n = new Noise("C:\\Users\\Василий\\Pictures\\2.jpg");
            n.Rayleigh(5, 10000).Save("C:\\Users\\Василий\\Pictures\\2.Noise.jpg", ImageFormat.Jpeg);
            Filters fil = new Filters();
            fil.Median(n.imageNew).Save("C:\\Users\\Василий\\Pictures\\2.Median.jpg", ImageFormat.Jpeg);
            Bitmap med = fil.resultImage;
            fil.Harmonic(n.imageNew).Save("C:\\Users\\Василий\\Pictures\\2.Harmonic.jpg", ImageFormat.Jpeg);
            Bitmap har = fil.resultImage;

            CompareImage first = new CompareImage(n.imageOld, med);
            Console.WriteLine(first.MSE(0, first.image2.Width, 0, first.image1.Height));
            Console.WriteLine(first.MSE_Part(10, 10));
            Console.WriteLine(first.YIK(0, first.image2.Width, 0, first.image1.Height));
            Console.WriteLine(first.YIK_Part(10, 10));

            Console.WriteLine("\n");

            CompareImage second = new CompareImage(n.imageOld, har);
            Console.WriteLine(second.MSE(0, second.image2.Width, 0, second.image1.Height));
            Console.WriteLine(second.MSE_Part(10, 10));
            Console.WriteLine(second.YIK(0, second.image2.Width, 0, second.image1.Height));
            Console.WriteLine(second.YIK_Part(10, 10));
        }
    }
    class CompareImage
    {
        public Bitmap image1 { get; set; }
        public Bitmap image2 { get; set; }
        public CompareImage(String im1, String im2) {
            image1 = new Bitmap(im1);
            image2 = new Bitmap(im2);
        }
        public CompareImage(Bitmap image1, Bitmap image2)
        {
            this.image1 = image1;
            this.image2 = image2;
        }

        public double MSE(int x1, int x2, int y1, int y2)
        {
            double res = 0;
            for(int i = x1; i < x2; i++)
                for (int j = y1; j < y2; j++)
                {
                    res += Math.Pow((image1.GetPixel(i, j).R - image2.GetPixel(i, j).R), 2);
                }
            return res / (x2 - x1) / (y2 - y1);
        }

        public double MSE_Part(int x, int y)
        {
            double res = 0;
            int k = 0, l = 0;
            for (k = 0; k < image1.Size.Width / x; k++)
                for (l = 0; l < image1.Size.Height / y; l++)
                {
                    res += MSE(k * x, (k + 1) * x, l * y, (l + 1) * y);
                }
            return res / l / k;
        }

        public double YIK(int x1, int x2, int y1, int y2)
        {
            double avgX = 0;
            double avgY = 0;
            double Qx = 0;
            double Qy = 0;
            double Qxy = 0;

            for (int i = x1; i < x2; i++)
                for (int j = y1; j < y2; j++)
                {
                    avgX += image1.GetPixel(i, j).R;
                    avgY += image2.GetPixel(i, j).R;
                }

            avgX /= (x2 - x1) * (y2 - y1);
            avgY /= (x2 - x1) * (y2 - y1);

            for (int i = x1; i < x2; i++)
                for (int j = y1; j < y2; j++)
                {
                    Qx += Math.Pow(image1.GetPixel(i, j).R - avgX, 2);
                    Qy += Math.Pow(image2.GetPixel(i, j).R - avgY, 2);
                    Qxy += (image1.GetPixel(i, j).R - avgX) * (image2.GetPixel(i, j).R - avgY);
                }
            Qx /= (x2 - x1) * (y2 - y1);
            Qxy /= (x2 - x1) * (y2 - y1);
            Qy /= (x2 - x1) * (y2 - y1);

            return (4 * Qxy + 0.001)*(avgX * avgY + 0.002) / (Qx + Qy + 0.001) / (Math.Pow(avgX, 2) + Math.Pow(avgY, 2) + 0.002);
        }

        public double YIK_Part(int x, int y)
        {
            double res = 0;
            int k = 0, l = 0;
            for (k = 0;  k < image1.Size.Width / x; k++)
                for (l = 0; l < image1.Size.Height / y; l++)
                {
                    res += YIK(k * x, (k + 1) * x, l * y, (l + 1) * y);
                }
            return res / l / k;
        }
    }

    class Noise
    {
        public Bitmap imageOld { get; set; }
        public Bitmap imageNew { get; }
        // шум райли (метод монте карло) и гармоническое среднее и медианый фильтр
        public Noise(String im) {
            imageOld = new Bitmap(im);
            imageNew = new Bitmap(im);
        }

        public Noise(Bitmap im)
        {
            this.imageOld = im;
            this.imageNew = im;
        }

        public Bitmap Rayleigh(double a, double b)
        {
            double[] p = new double[256]; // распределение вероятностей интенсивности шума на отрезок [0,1]
            for (int z = 0; z < 256; z++ ) // заполнение p по формуле
            {
                if (z < a)
                    p[z] = 0;
                else
                    p[z] = p[z - 1] + 2/b * (z - a) * Math.Exp(-(z - a) * (z - a) / b);
            }
            //метод Монте-Карло 
            Random random = new Random();
            for (int i = 0; i < imageOld.Width; i++)
                for (int j = 0; j < imageOld.Height; j++)
                {
                    if(random.Next(100) <= 5)
                    {
                        double noise = random.NextDouble();
                        for (int k = 0; k < 256; k++)
                        {
                            if (noise < p[k])
                            {
                                imageNew.SetPixel(i, j, Color.FromArgb(k - 1, k - 1, k - 1));
                                break;
                            }
                        }
                    }
                    else
                        imageNew.SetPixel(i, j, imageOld.GetPixel(i, j));
                }
            return imageNew;
        }

        
    }

    class Filters
    {
        public Bitmap resultImage;
        public Bitmap Median(string Name)
        {
            Bitmap sourceImage = new Bitmap(Name);
            resultImage = new Bitmap(sourceImage.Width, sourceImage.Height);
            for (int i = 0; i < sourceImage.Width; i++)
                for (int j = 0; j < sourceImage.Height; j++)
                    resultImage.SetPixel(i, j, calculateNewPixelColorMedian(sourceImage, i, j));
            return resultImage;
        }
        public Bitmap Median(Bitmap im)
        {
            Bitmap sourceImage = im;
            resultImage = new Bitmap(sourceImage.Width, sourceImage.Height);
            for (int i = 0; i < sourceImage.Width; i++)
                for (int j = 0; j < sourceImage.Height; j++)
                    resultImage.SetPixel(i, j, calculateNewPixelColorMedian(sourceImage, i, j));
            return resultImage;
        }
        protected Color calculateNewPixelColorMedian(Bitmap sourceImage, int i, int j)
        {
            int radiusX = 3 / 2;
            int radiusY = 3 / 2;
            int q = 0;
            int[] colorsR = new int[9];
            int[] colorsG = new int[9];
            int[] colorsB = new int[9];
            for (int l = -radiusY; l <= radiusY; l++)
                for (int k = -radiusX; k <= radiusX; k++)
                {
                    int idX = Clamp(i + k, 0, sourceImage.Width - 1);
                    int idY = Clamp(j + l, 0, sourceImage.Height - 1);
                    Color neighborColor = sourceImage.GetPixel(idX, idY);
                    colorsR[q] = neighborColor.R;
                    colorsG[q] = neighborColor.G;
                    colorsB[q] = neighborColor.B;
                    q++;
                }
            Sort(colorsR);
            Sort(colorsB);
            Sort(colorsG);
            return Color.FromArgb(colorsR[9 / 2], colorsG[9 / 2], colorsB[9 / 2]);
        }

        public Bitmap Harmonic(string Name)
        {
            Bitmap sourceImage = new Bitmap(Name);
            resultImage = new Bitmap(sourceImage.Width, sourceImage.Height);
            for (int i = 0; i < sourceImage.Width; i++)
                for (int j = 0; j < sourceImage.Height; j++)
                    resultImage.SetPixel(i, j, calculateNewPixelColorHarmonic(sourceImage, i, j));
            return resultImage;
        }
        public Bitmap Harmonic(Bitmap im)
        {
            Bitmap sourceImage = im;
            resultImage = new Bitmap(sourceImage.Width, sourceImage.Height);
            for (int i = 0; i < sourceImage.Width; i++)
                for (int j = 0; j < sourceImage.Height; j++)
                    resultImage.SetPixel(i, j, calculateNewPixelColorHarmonic(sourceImage, i, j));
            return resultImage;
        }
        protected Color calculateNewPixelColorHarmonic(Bitmap sourceImage, int i, int j)
        {
            int radiusX = 3 / 2;
            int radiusY = 3 / 2;
            double colorsR = 0;
            double colorsG = 0;
            double colorsB = 0;
            for (int l = -radiusY; l <= radiusY; l++)
                for (int k = -radiusX; k <= radiusX; k++)
                {
                    int idX = Clamp(i + k, 0, sourceImage.Width - 1);
                    int idY = Clamp(j + l, 0, sourceImage.Height - 1);
                    Color neighborColor = sourceImage.GetPixel(idX, idY);
                    if(neighborColor.R != 0)
                        colorsR += 1f / neighborColor.R;
                    if (neighborColor.G != 0)
                        colorsG += 1f  /neighborColor.G;
                    if (neighborColor.B != 0)
                        colorsB += 1f / neighborColor.B;
                }
            return Color.FromArgb(Clamp((int)(9 / colorsR), 0, 255), Clamp((int)(9 / colorsG), 0, 255), Clamp((int)(9 / colorsB), 0, 255));
        }
        private void Swap(ref int x, ref int y)
        {
            int t = x;
            x = y;
            y = t;
        }
        private void Sort(int[] arr)
        {
            for (int i = 0; i < arr.Length; i++)
                for (int j = 0; j < arr.Length; j++)
                {
                    if (arr[i] > arr[j])
                        Swap(ref arr[i], ref arr[j]);
                }
        }
        public int Clamp(int value, int min, int max)
        {
            if (value > max) return max;
            if (value < min) return min;
            return value;
        }
    }
}
