// See https://aka.ms/new-console-template for more information

using System;
using System.ComponentModel;
using System.ComponentModel.Design;
using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.CompilerServices;
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
            

            Noise n = new Noise("C:\\Users\\Василий\\Pictures\\2.jpg");
            n.Rayleigh(5, 10000).Save("C:\\Users\\Василий\\Pictures\\2.Noise.jpg", ImageFormat.Jpeg);
            Filters fil = new Filters();
            fil.Median(n.imageNew).Save("C:\\Users\\Василий\\Pictures\\2.Median.jpg", ImageFormat.Jpeg);
            Bitmap med = fil.resultImage;
            fil.Harmonic(n.imageNew).Save("C:\\Users\\Василий\\Pictures\\2.Harmonic.jpg", ImageFormat.Jpeg);
            Bitmap har = fil.resultImage;

            CompareImage noise = new CompareImage(n.imageOld, n.imageNew);
            Console.WriteLine(noise.MSE(0, noise.image2.Width, 0, noise.image1.Height));
            Console.WriteLine(noise.MSE_Part(10, 10));
            Console.WriteLine(noise.YIK(0, noise.image2.Width, 0, noise.image1.Height));
            Console.WriteLine(noise.YIK_Part(10, 10));

            Console.WriteLine("\n");

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
            */

            Bitmap image = new Bitmap("C:\\Users\\Василий\\Pictures\\5.1.png");
            Borders.Kenny(image).Save("C:\\Users\\Василий\\Pictures\\5.1.Kenny.png", ImageFormat.Png);
            Bitmap image2 = new Bitmap("C:\\Users\\Василий\\Pictures\\6.jpg");
            //Borders.Moments(image2);

            /*Bitmap res = WaveletTransform.HaarWaveletTransform(image);
            res.Save("C:\\Users\\Василий\\Pictures\\5.Wavelet.png", ImageFormat.Png);
            WaveletTransform.InverseHaarWaveletTransform(res).Save("C:\\Users\\Василий\\Pictures\\5.ReserveWavelet.png", ImageFormat.Png);
            WaveletTransform.LowPassFilter(image).Save("C:\\Users\\Василий\\Pictures\\5.LowPass.png", ImageFormat.Png);
            WaveletTransform.HighPassFilter(image).Save("C:\\Users\\Василий\\Pictures\\5.HighPass.png", ImageFormat.Png);*/
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
                {
                    p[z] = 0;
                    //Console.WriteLine(p[z]);
                }
                else
                {
                    p[z] = p[z - 1] + 2 / b * (z - a) * Math.Exp(-(z - a) * (z - a) / b);
                   //Console.WriteLine(2 / b * (z - a) * Math.Exp(-(z - a) * (z - a) / b));
                }
            }
            Console.WriteLine("\n");
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

    static class Borders
    {
        static public Bitmap Kenny(Bitmap imageOld)
        {
            Bitmap smoothImage = GaussianFilter(imageOld, 3, 3); // сглаживание изображения (удаление шума Гауссом)
            double[,] angles;
            double[,] sobelImage = Sobel(smoothImage, out angles); // вычистление значения и направления градиентов
            double[,] maxImage = NoMax(sobelImage, ref angles); // подавление немаксимумов (оставляем только локальные максимумы в направлении градиента)
            Bitmap thresholdIm = Threshold(maxImage, 0.05f, 0.55f); // подавление несуществующих контуров (сравниваем с пороговыми значениями)
            Bitmap ambiguityIm = Ambiguity(thresholdIm); // смотриим связанность оставшихся пикселей (если они они не соприкасаются с отсальными пикселями - удаляем)
            thresholdIm.Save("C:\\Users\\Василий\\Pictures\\5.1.Екуыр.png", ImageFormat.Png);
            return ambiguityIm;
        }
        //толстая граница и что-то с моментами дллжно быть много
        static private Bitmap Ambiguity(Bitmap image) {
            Bitmap newIm = new Bitmap(image);
            int changes = 0;
            do
            {
                changes = 0;
                for (int x = 0; x < image.Width; x++)
                    for (int y = 0; y < image.Height; y++)
                    {
                        if (newIm.GetPixel(x, y).R == 127)
                        {
                            bool flag = false;
                            for (int dx = -1; dx <= 1; dx++)
                            {
                                for (int dy = -1; dy <= 1; dy++)
                                {
                                    int newX = x + dx;
                                    int newY = y + dy;
                                    if (newX >= 0 && newX < image.Width && newY >= 0 && newY < image.Height)
                                    {
                                        Color neighborColor = newIm.GetPixel(newX, newY);

                                        if (neighborColor.R > 127)
                                        {
                                            newIm.SetPixel(x, y, Color.FromArgb(255, 255, 255));
                                            changes++;
                                            flag = true;
                                            break;
                                        }
                                    }
                                }
                                if (flag)
                                    break;

                            }
                            /*if (!flag)
                            {
                                newIm.SetPixel(x, y, Color.FromArgb(0, 0, 0));
                            }*/
                        }
                    }
            } while (changes != 0);
            for (int x = 0; x < image.Width; x++)
                for (int y = 0; y < image.Height; y++)
                {
                    if (newIm.GetPixel(x, y).R == 127)
                    {
                        newIm.SetPixel(x, y, Color.FromArgb(0, 0, 0));
                    }
                }
           return newIm;
        }
        static private Bitmap Threshold(double[,] image, float low_p, float high_p) { // пороги в долях от 0 до 1
            double max = 0.0;
            for (int x = 0; x < image.GetLength(0); x++)
                for (int y = 0; y < image.GetLength(1); y++)
                {
                    if (image[x, y] > max)
                        max = image[x, y];
                }
            int low = (int)(low_p * max);
            int high = (int)(high_p * max);
            
            Bitmap newIm = new Bitmap(image.GetLength(0), image.GetLength(1));

            for (int x = 0; x < image.GetLength(0); x++)
                for (int y = 0; y < image.GetLength(1); y++)
                {
                    if (image[x, y] > high)
                        newIm.SetPixel(x, y, Color.FromArgb(255, 255, 255));
                    else if (image[x, y] < low)
                        newIm.SetPixel(x, y, Color.FromArgb(0, 0, 0));
                    else
                        newIm.SetPixel(x, y, Color.FromArgb(127, 127, 127));
                }
            return newIm;
        }
        static double[,] NoMax(double[,] image, ref double[,] angles) // тут что-то не так
        {
            //Bitmap newIm = new Bitmap(image); //356 397
            double[,] res = (double[,])image.Clone();
            for (int x = 0; x < image.GetLength(0); x++)
            {
                for (int y = 0; y < image.GetLength(1); y++)
                {
                    double angle = angles[x, y];

                    int neighborX1 = 0, neighborY1 = 0, neighborX2 = 0, neighborY2 = 0;

                    if ((angle >= - Math.PI / 8 && angle < Math.PI / 8) || (angle >= 7 * Math.PI / 8 && angle <= 8 * Math.PI / 8) || (angle >= -8 * Math.PI / 8 && angle <= -7 * Math.PI / 8))
                    {
                        neighborX1 = x;
                        neighborY1 = y + 1;
                        neighborX2 = x;
                        neighborY2 = y - 1;
                    }
                    else if ((angle >= Math.PI / 8 && angle < 3 * Math.PI / 8) || (angle >= -7 * Math.PI / 8 && angle < -5 * Math.PI / 8))
                    {
                        neighborX1 = x - 1;
                        neighborY1 = y - 1;
                        neighborX2 = x + 1;
                        neighborY2 = y + 1;
                    }
                    else if ((angle >= 3 * Math.PI / 8 && angle < 5 * Math.PI / 8) || (angle >= -5 * Math.PI / 8 && angle < -3 * Math.PI / 8))
                    {
                        neighborX1 = x + 1;
                        neighborY1 = y;
                        neighborX2 = x - 1;
                        neighborY2 = y;
                    }
                    else if ((angle >= 5 * Math.PI / 8 && angle < 7 * Math.PI / 8) || (angle >= -3 * Math.PI / 8 && angle < -1 * Math.PI / 8))
                    {
                        neighborX1 = x + 1;
                        neighborY1 = y - 1;
                        neighborX2 = x - 1;
                        neighborY2 = y + 1;
                    }
                    if (neighborX1 < 0 || neighborX1 >= image.GetLength(0) || neighborY1 < 0 || neighborY1 >= image.GetLength(1) ||
                        neighborX2 < 0 || neighborX2 >= image.GetLength(0) || neighborY2 < 0 || neighborY2 >= image.GetLength(1))
                    {
                        continue;
                    }

                    if (image[x, y] <= image[neighborX1, neighborY1] || image[x, y] <= image[neighborX2, neighborY2])
                    {
                        res[x, y] = 0;
                    }

                }
            }
            return res;

        }
        static private double[,] Sobel(Bitmap sourceImage, out double[,] angles)
        {
            int[,] kernelX = {
            { -1, 0, 1 },
            { -2, 0, 2 },
            { -1, 0, 1 }
        };

            int[,] kernelY = {
            { -1, -2, -1 },
            { 0, 0, 0 },
            { 1, 2, 1 }
        };
            
            int width = sourceImage.Width;
            int height = sourceImage.Height;
            double[,] res = new double[width, height];
            Bitmap outputImage = new Bitmap(width, height);
            angles = new double[width, height];

            for (int x = 1; x < width - 1; x++)
            {
                for (int y = 1; y < height - 1; y++)
                {
                    int gx = 0, gy = 0;

                    for (int i = -1; i <= 1; i++)
                    {
                        for (int j = -1; j <= 1; j++)
                        {
                            int idX1 = Clamp(x + i, 0, sourceImage.Width - 1);
                            int idY1 = Clamp(y + j, 0, sourceImage.Height - 1);

                            Color pixel = sourceImage.GetPixel(idX1, idY1);
                            int gray = pixel.R;
                            //int gray = (int)(pixel.R * 0.299 + pixel.G * 0.587 + pixel.B * 0.114);

                            gx += gray * kernelX[i + 1, j + 1];
                            gy += gray * kernelY[i + 1, j + 1];
                        }
                    }

                    double magnitude = Math.Sqrt(gx * gx + gy * gy);

                    //magnitude = Math.Min(255, Math.Max(0, magnitude));
                    res[x, y] = magnitude;
                    double angle = Math.Atan2(gy, gx);

                    angles[x, y] = angle;

                    //outputImage.SetPixel(x, y, Color.FromArgb(magnitude, magnitude, magnitude));
                }
            }

            return res;
        }
        static private int Clamp(int value, int min, int max)
        {
            if (value > max) return max;
            if (value < min) return min;
            return value;
        }
        static private Bitmap GaussianFilter(Bitmap image, int radius, float sigma)
        {
            int size = 2 * radius + 1;
            float[,] kernel = new float[size, size];
            float norm = 0;
            for (int i = -radius; i <= radius; i++)
                for (int j = -radius; j <= radius; j++)
                {
                    kernel[i + radius, j + radius] = (float)(Math.Exp(-(i * i + j * j) / (sigma * sigma)));
                    norm += kernel[i + radius, j + radius];
                }
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    kernel[i, j] /= norm;

            Bitmap resultImage = new Bitmap(image.Width, image.Height);
            var g = Graphics.FromImage(resultImage);
            g.Clear(Color.White);
            for (int i = 0; i < image.Width; i++)
            {
                for (int j = 0; j < image.Height; j++)
                    resultImage.SetPixel(i, j, NewPixelGaussian(image, i, j, ref kernel));
            }
            return resultImage;
        }
        static private Color NewPixelGaussian(Bitmap sourceImage, int i, int j, ref float[,] kernel)
        {
            int radiusX = kernel.GetLength(0) / 2;
            int radiusY = kernel.GetLength(1) / 2;
            float resultR = 0;
            for (int l = -radiusY; l <= radiusY; l++)
                for (int k = -radiusX; k <= radiusX; k++)
                {
                    int idX = Clamp(i + k, 0, sourceImage.Width - 1);
                    int idY = Clamp(j + l, 0, sourceImage.Height - 1);
                    Color neighborColor = sourceImage.GetPixel(idX, idY);
                    resultR += neighborColor.R * kernel[k + radiusX, l + radiusY];
                }
            return Color.FromArgb(Clamp((int)resultR, 0, 255), Clamp((int)resultR, 0, 255), Clamp((int)resultR, 0, 255));
        }


        static public void Moments(Bitmap im)
        {
            Bitmap image = new Bitmap(im);

            int[,] objs = new int[image.Width, image.Height];
            int objNum = 1;
            Array.Clear(objs);

            Bitmap binaryImage = ApplyThreshold(im, 128);

            for (int x = 0; x < binaryImage.Width; x++)
            {
                for (int y = 0; y < binaryImage.Height; y++)
                {
                    if (binaryImage.GetPixel(x, y).R == 255 && objs[x, y] == 0)
                    {
                        int a = binaryImage.GetPixel(x, y).A;
                        SearchObj(binaryImage, objs, x, y, objNum);
                        objNum++;
                    }
                }
            }

            for (int i = 1; i < objNum; i++)
            {
                Console.WriteLine("Object " + i);
                Console.WriteLine("Simple moments");
                for (int p = 0; p < 3; p++)
                    for (int q = 0; q < 3; q++)
                    {
                        Console.WriteLine("M" + p + q + " = " + FindSimpleMoment(objs, i, p, q));
                    }
                Console.WriteLine("Central moments");
                for (int p = 0; p < 3; p++)
                    for (int q = 0; q < 3; q++)
                    {
                        Console.WriteLine("M" + p + q + " = " + FindCentralMoment(objs, i, p, q));
                    }
            }

        }
        static private void SearchObj(Bitmap image, int[,] arr, int sx, int sy, int objNum)
        {
            arr[sx, sy] = objNum;
            List<List<int>> saves = new List<List<int>>
            {
                new List<int> { sx, sy }
            };
            while (saves.Count != 0)
            {
                List<int> p = saves[0];
                int x = p[0];
                int y = p[1];
                for (int i = -1; i < 2; i++)
                {
                    for (int j = -1; j < 2; j++)
                    {
                        if (x + i >= arr.GetLength(0) || y + j >= arr.GetLength(1) || x + i < 0 || y + j < 0)
                            continue;
                        if (arr[x + i, y + j] == 0 && image.GetPixel(x + i, y + j).R == 255)
                        {
                            arr[x + i, y + j] = objNum;
                            saves.Add(new List<int> { x + i, y + j });
                            //SearchObj(image, arr, x + i, y + j, objNum);
                        }
                    }
                }
                saves.Remove(p);
            }
        }
        static private Bitmap ApplyThreshold(Bitmap grayImage, int threshold)
        {
            Bitmap binaryImage = new Bitmap(grayImage.Width, grayImage.Height);

            for (int x = 0; x < grayImage.Width; x++)
            {
                for (int y = 0; y < grayImage.Height; y++)
                {
                    Color pixel = grayImage.GetPixel(x, y);
                    int value = pixel.R;
                    Color binaryColor = (value < threshold) ? Color.Black : Color.White;
                    binaryImage.SetPixel(x, y, binaryColor);
                }
            }
            return binaryImage;
        }
        static private double FindSimpleMoment(int[,] objs, int objNum, int p, int q)
        {
            double res = 0;
            for (int i = 0; i < objs.GetLength(0); i++)
                for (int j = 0; j < objs.GetLength(1); j++)
                {
                    if (objs[i, j] == objNum)
                    {
                        res += Math.Pow(i, p) * Math.Pow(j, q);
                    }
                }
            return res;
        }
        static private double FindCentralMoment(int[,] objs, int objNum, int p, int q)
        {
            double res = 0;
            double M00 = FindSimpleMoment(objs, objNum, 0, 0);
            double M10 = FindSimpleMoment(objs, objNum, 1, 0);
            double M01 = FindSimpleMoment(objs, objNum, 0, 1);

            double xCent = M10 / M00;
            double yCent = M01 / M00;

            for (int i = 0; i < objs.GetLength(0); i++)
                for (int j = 0; j < objs.GetLength(1); j++)
                {
                    if (objs[i, j] == objNum)
                    {
                        res += Math.Pow(i - xCent, p) * Math.Pow(j - yCent, q);
                    }
                }
            return res;
        }
    }

    static class WaveletTransform
    {
        static double[,] lastArr;
        static public Bitmap HaarWaveletTransform(Bitmap im)
        {
            double[,] image = LoadImageToData(im);
            int rows = image.GetLength(0);
            int cols = image.GetLength(1);

            for (int i = 0; i < rows; i++)
            {
                double[] row = new double[cols];
                for (int j = 0; j < cols; j++)
                {
                    row[j] = image[i, j];
                }
                HaarWaveletTransform1D(row);

                for (int j = 0; j < cols; j++)
                {
                    image[i, j] = row[j];
                }
            }

            for (int j = 0; j < cols; j++)
            {
                double[] column = new double[rows];
                for (int i = 0; i < rows; i++)
                {
                    column[i] = image[i, j];
                }
                HaarWaveletTransform1D(column);

                for (int i = 0; i < rows; i++)
                {
                    image[i, j] = column[i];
                }
            }
            lastArr = (double[,])image.Clone();
            return SaveDataToImage(image);
        }
        static void HaarWaveletTransform1D(double[] data)
        {
            int length = data.Length;

            int half = length / 2;
            double[] temp = new double[length];

            for (int i = 0, j = 0; i < half; i++, j += 2)
            {
                temp[i] = (data[j] + data[j + 1]) / 2;
                temp[i + half] = (data[j] - data[j + 1]) / 2;
            }

            Array.Copy(temp, data, length);
        }
        static public Bitmap InverseHaarWaveletTransform(Bitmap im)
        {
            double[,] image = lastArr;
            int rows = image.GetLength(0);
            int cols = image.GetLength(1);

            for (int j = 0; j < cols; j++)
            {
                double[] column = new double[rows];
                for (int i = 0; i < rows; i++)
                {
                    column[i] = image[i, j];
                }
                InverseHaarWaveletTransform1D(column);

                for (int i = 0; i < rows; i++)
                {
                    image[i, j] = column[i];
                }
            }

            for (int i = 0; i < rows; i++)
            {
                double[] row = new double[cols];
                for (int j = 0; j < cols; j++)
                {
                    row[j] = image[i, j];
                }
                InverseHaarWaveletTransform1D(row);

                for (int j = 0; j < cols; j++)
                {
                    image[i, j] = row[j];
                }
            }
            return SaveDataToImageInverse(image);
        }

        static void InverseHaarWaveletTransform1D(double[] data)
        {
            int length = data.Length;
            int half = length / 2;

            double[] temp = new double[length];

            for (int i = 0, j = 0; i < half; i++, j += 2)
            {
                temp[j] = (data[i] + data[i + half]) ;
                temp[j + 1] = (data[i] - data[i + half]) ;
            }

            Array.Copy(temp, data, length);
        }

        static public Bitmap LowPassFilter(Bitmap im)
        {
            double[,] image = LoadImageToData(im);
            int rows = image.GetLength(0);
            int cols = image.GetLength(1);

            for (int i = 0; i < rows; i++)
            {
                double[] row = new double[cols];
                for (int j = 0; j < cols; j++)
                {
                    row[j] = image[i, j];
                }
                LowPassFilter1D(row);

                for (int j = 0; j < cols; j++)
                {
                    image[i, j] = row[j];
                }
            }

            for (int j = 0; j < cols; j++)
            {
                double[] column = new double[rows];
                for (int i = 0; i < rows; i++)
                {
                    column[i] = image[i, j];
                }
                LowPassFilter1D(column);

                for (int i = 0; i < rows; i++)
                {
                    image[i, j] = column[i];
                }
            }

            return SaveDataToImage(image);
        }
        static void LowPassFilter1D(double[] data)
        {
            int length = data.Length;
            int half = length / 2;

            double[] temp = new double[length];

            for (int i = 0, j = 0; i < half; i++, j += 2)
            {
                temp[j] = (data[j] + data[j + 1]) / 2;
                temp[j + 1] = (data[j] + data[j + 1]) / 2;
            }

            Array.Copy(temp, data, length);
        }
        static public Bitmap HighPassFilter(Bitmap im)
        {
            double[,] image = LoadImageToData(im);
            int rows = image.GetLength(0);
            int cols = image.GetLength(1);

            for (int i = 0; i < rows; i++)
            {
                double[] row = new double[cols];
                for (int j = 0; j < cols; j++)
                {
                    row[j] = image[i, j];
                }
                HighPassFilter1D(row);

                for (int j = 0; j < cols; j++)
                {
                    image[i, j] = row[j];
                }
            }

            for (int j = 0; j < cols; j++)
            {
                double[] column = new double[rows];
                for (int i = 0; i < rows; i++)
                {
                    column[i] = image[i, j];
                }
                HighPassFilter1D(column);

                for (int i = 0; i < rows; i++)
                {
                    image[i, j] = column[i];
                }
            }

            return SaveDataToImage(image);
        }
        static void HighPassFilter1D(double[] data)
        {
            int length = data.Length;
            int half = length / 2;

            double[] temp = new double[length];

            for (int i = 0, j = 0; i < half; i++, j += 2)
            {
                temp[j] = (data[j] - data[j + 1]) / 2;
                temp[j + 1] = (data[j] - data[j + 1]) / 2;
            }

            Array.Copy(temp, data, length);
        }
        static double[,] LoadImageToData(Bitmap image)
        {
            int width = image.Width;
            int height = image.Height;

            double[,] data = new double[width, height];

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    Color pixelColor = image.GetPixel(x, y);
                    data[x, y] = pixelColor.R;
                }
            }

            return data;
        }
        static double[,] LoadImageToDataInverse(Bitmap image)
        {
            int width = image.Width;
            int height = image.Height;

            double[,] data = new double[width, height];

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    int pixelColor = image.GetPixel(x, y).R;
                    if (x < width / 2)
                    {
                        if (y < height / 2)
                        {
                            data[x, y] = pixelColor;
                        }
                        else
                        {
                            data[x, y] = pixelColor - 128;
                        }
                    }
                    else
                    {
                        data[x, y] = pixelColor - 128;
                    }
                }
            }

            return data;
        }

        static Bitmap SaveDataToImage(double[,] data)
        {
            int width = data.GetLength(0);
            int height = data.GetLength(1);

            Bitmap resultImage = new Bitmap(width, height);

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    int pixelValue;
                    if (data[x, y] < 0)
                        pixelValue = (int)data[x, y] + 128;
                    else
                        pixelValue = (int)data[x, y];
                    double a = data[x, y];
                    Color pixelColor = Color.FromArgb(pixelValue, pixelValue, pixelValue);
                    resultImage.SetPixel(x, y, pixelColor);
                }
            }

            return resultImage;
        }
        static Bitmap SaveDataToImageInverse(double[,] data)
        {
            int width = data.GetLength(0);
            int height = data.GetLength(1);

            Bitmap resultImage = new Bitmap(width, height);

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    int pixelValue = (int)(data[x, y]);
                    Color pixelColor = Color.FromArgb(pixelValue, pixelValue, pixelValue);
                    resultImage.SetPixel(x, y, pixelColor);
                }
            }

            return resultImage;
        }
    }
}
// Выделение границ и поиск объектов при помощи методов Хаффа
// метод кени и моменты