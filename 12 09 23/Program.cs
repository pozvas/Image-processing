// See https://aka.ms/new-console-template for more information

using System.ComponentModel;
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

            Bitmap image = new Bitmap("C:\\Users\\Василий\\Pictures\\5.Cenny.png");
            //Borders.Cenny(image).Save("C:\\Users\\Василий\\Pictures\\5.Cenny.png", ImageFormat.Png);
            Borders.Haff(image).Save("C:\\Users\\Василий\\Pictures\\5.Haff.png", ImageFormat.Png);
            
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
        static private float[,] kernelY = { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } };
        static private float[,] kernelX = { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } };
        static public Bitmap Cenny(Bitmap imageOld)
        {
            Bitmap smoothImage = GaussianFilter(imageOld, 7, 1); // сглаживание изображения (удаление шума Гауссом)
            double[,] angles = new double[imageOld.Width, imageOld.Height];
            Bitmap sobelImage = Sobel(smoothImage, ref angles); // вычистление значения и направления градиентов
            Bitmap maxImage = NoMax(sobelImage, ref angles); // подавление немаксимумов (оставляем только локальные максимумы в направлении градиента)
            Bitmap thresholdIm = Threshold(maxImage, 0.45f, 0.6f); // подавление несуществующих контуров (сравниваем с пороговыми значениями)
            Bitmap ambiguityIm = Ambiguity(thresholdIm); // смотриим связанность оставшихся пикселей (если они они не соприкасаются с отсальными пикселями - удаляем)
            return ambiguityIm;
        }
        static private Bitmap Ambiguity(Bitmap image) {
            Bitmap newIm = new Bitmap(image);
            for (int x = 0; x < image.Width; x++)
                for (int y = 0; y < image.Height; y++)
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
                                Color neighborColor = image.GetPixel(newX, newY);

                                // Если соседний пиксель также является граничным и его интенсивность превышает порог,
                                // свяжите оба пикселя в контур
                                if (neighborColor.R > 127)
                                {
                                    newIm.SetPixel(x, y, Color.FromArgb(255, 255, 255));
                                    flag = true;
                                    break;
                                }
                            }
                        }
                        if (flag)
                            break;
                        else
                        {
                            newIm.SetPixel(x, y, Color.FromArgb(0, 0, 0));
                        }
                    }
                }
            return newIm;
        }
        static private Bitmap Threshold(Bitmap image, float low_p, float high_p) { // пороги в долях от 0 до 1
            int low = (int)(low_p * 255f);
            int high = (int)(high_p * 255f);
            Bitmap newIm = new Bitmap(image.Width, image.Height);
            for (int x = 0; x < image.Width; x++)
                for (int y = 0; y < image.Height; y++)
                {
                    if (image.GetPixel(x, y).B > high)
                        newIm.SetPixel(x, y, Color.FromArgb(255, 255, 255));
                    else if (image.GetPixel(x, y).B < low)
                        newIm.SetPixel(x, y, Color.FromArgb(0, 0, 0));
                    else
                        newIm.SetPixel(x, y, Color.FromArgb(127, 127, 127));
                }
            return newIm;
        }
        static private Bitmap NoMax(Bitmap image, ref double[,] angles)
        {
            Bitmap newIm = new Bitmap(image); //356 397
            for (int x = 0; x < image.Width; x++)
                for (int y = 0; y < image.Height; y++)
                {
                    int neighborX1 = (int)Math.Round(x + Math.Cos(angles[x, y]));
                    int neighborY1 = (int)Math.Round(y + Math.Sin(angles[x, y]));
                    int neighborX2 = (int)Math.Round(x - Math.Cos(angles[x, y]));
                    int neighborY2 = (int)Math.Round(y - Math.Sin(angles[x, y]));

                    bool first = !(neighborY1 < 0 || neighborY1 >= image.Height || neighborX1 < 0 || neighborX1 >= image.Width);
                    bool second = !(neighborY2 < 0 || neighborY2 >= image.Height || neighborX2 < 0 || neighborX2 >= image.Width);

                    if (first && second)
                    {
                        if (!(image.GetPixel(x, y).B >= image.GetPixel(neighborX1, neighborY1).B && image.GetPixel(x, y).B >= image.GetPixel(neighborX2, neighborY2).B))
                        {
                            newIm.SetPixel(x, y, Color.Black);
                        }
                    }
                    else if (!first && second)
                    {
                        if (!(image.GetPixel(x, y).B >= image.GetPixel(neighborX2, neighborY2).B))
                        {
                            newIm.SetPixel(x, y, Color.Black);
                        }
                    }
                    else if (first && !second)
                    {
                        if (!(image.GetPixel(x, y).B >= image.GetPixel(neighborX1, neighborY1).B))
                        {
                            newIm.SetPixel(x, y, Color.Black);
                        }
                    }
                    else
                    {
                        continue;
                    }

                }
            return newIm;
                    
        }
        static private Bitmap Sobel(Bitmap sourceImage, ref double[,] angles)
        {
            Bitmap resultImage = new Bitmap(sourceImage.Width, sourceImage.Height);
            var g = Graphics.FromImage(resultImage);
            g.Clear(Color.White);
            for (int i = 0; i < sourceImage.Width; i++)
            {
                for (int j = 0; j < sourceImage.Height; j++)
                    resultImage.SetPixel(i, j, NewPixelSobel(sourceImage, i, j, ref angles[i,j]));
            }
            return resultImage;
        }
        static private int Clamp(int value, int min, int max)
        {
            if (value > max) return max;
            if (value < min) return min;
            return value;
        }
        static private Color NewPixelSobel(Bitmap sourceImage, int i, int j, ref double angle)
        {
            float resultRY = 0, resultGY = 0, resultBY = 0;
            int radiusX1 = kernelY.GetLength(0) / 2;
            int radiusY1 = kernelY.GetLength(1) / 2;

            float resultRX = 0, resultGX = 0, resultBX = 0;
            int radiusX2 = kernelX.GetLength(0) / 2;
            int radiusY2 = kernelX.GetLength(1) / 2;


            for (int l = -radiusY1; l <= radiusY1; l++)
            {
                for (int k = -radiusX1; k <= radiusX1; k++)
                {
                    int idX1 = Clamp(i + k, 0, sourceImage.Width - 1);
                    int idY1 = Clamp(j + l, 0, sourceImage.Height - 1);

                    Color neighborColor = sourceImage.GetPixel(idX1, idY1);

                    resultRY += neighborColor.R * kernelY[k + radiusX1, l + radiusY1];
                    //resultGY += neighborColor.G * kernelY[k + radiusX1, l + radiusY1];
                    //resultBY += neighborColor.B * kernelY[k + radiusX1, l + radiusY1];
                }
            }

            for (int l = -radiusY2; l <= radiusY2; l++)
            {
                for (int k = -radiusX2; k <= radiusX2; k++)
                {
                    int idX2 = Clamp(i + k, 0, sourceImage.Width - 1);
                    int idY2 = Clamp(j + l, 0, sourceImage.Height - 1);

                    Color neighborColor = sourceImage.GetPixel(idX2, idY2);

                    resultRX += neighborColor.R * kernelX[k + radiusX2, l + radiusY2];
                    //resultGX += neighborColor.G * kernelX[k + radiusX2, l + radiusY2];
                   // resultBX += neighborColor.B * kernelX[k + radiusX2, l + radiusY2];
                }
            }
            float resultR = (float)Math.Sqrt(Math.Pow(resultRX, 2) + Math.Pow(resultRY, 2));
            //float resultG = (float)Math.Sqrt(Math.Pow(resultGX, 2) + Math.Pow(resultGY, 2));
            //float resultB = (float)Math.Sqrt(Math.Pow(resultGX, 2) + Math.Pow(resultGY, 2));

            angle = Math.Round(Math.Atan2(resultRX, resultRY) / Math.PI * 4d) * Math.PI / 4; 
            //angle = Math.Atan2(resultBX, resultBY);

            return Color.FromArgb(Clamp((int)resultR, 0, 255), Clamp((int)resultR, 0, 255), Clamp((int)resultR, 0, 255));
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


        static public Bitmap Haff(Bitmap im)
        {
            Bitmap image = new Bitmap(im);

            Bitmap binaryImage = ApplyThreshold(im, 128);

            Point centerOfMass = CalculateCenterOfMass(binaryImage);

            using (Graphics g = Graphics.FromImage(image))
            {
                g.DrawEllipse(Pens.Red, centerOfMass.X - 5, centerOfMass.Y - 5, 10, 10);
                image.SetPixel(centerOfMass.X, centerOfMass.Y, Color.Red);
            }

            return image;
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
        static private Point CalculateCenterOfMass(Bitmap binaryImage)
        {
            int sumX = 0;
            int sumY = 0;
            int count = 0;

            for (int x = 0; x < binaryImage.Width; x++)
            {
                for (int y = 0; y < binaryImage.Height; y++)
                {
                    Color pixel = binaryImage.GetPixel(x, y);
                    if (pixel.R == 255) 
                    {
                        sumX += x;
                        sumY += y;
                        count++;
                    }
                }
            }

            if (count > 0)
            {
                int centerX = sumX / count;
                int centerY = sumY / count;
                return new Point(centerX, centerY);
            }
            else
            {
                return Point.Empty; 
            }
        }
    }
}
// Выделение границ и поиск объектов при помощи методов Хаффа
// метод кени и моменты