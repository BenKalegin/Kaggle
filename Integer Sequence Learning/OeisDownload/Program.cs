using System;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Threading;
using System.Threading.Tasks;

namespace OeisDownload
{
    public class Sequence
    {
        public string Id { get; set; }
        public string Seq { get; set; }
    }

    class Program
    {
        private const string inputTestCsv = @"../../../input/test.csv";
        private const string downloadFolder = @"../../../download";

        static string MakeUrl(string seq)
        {
            return "http://oeis.org/search?fmt=json&q=" + seq;

        }

        private static void DownloadFiles(Sequence[] items)
        {
            var todo = items.ToDictionary(i => i.Id, i => i);

            while (todo.Any())
            {
                var batch = todo.Take(40).Select(d => d.Value).ToList();
                int i = 0;
                while(i < batch.Count)
                {
                    var item = batch[i];
                    if (File.Exists(FileName(item)))
                    {
                        todo.Remove(item.Id);
                        batch.Remove(item);
                    }
                    else
                    {
                        i++;
                    }
                }

                if (batch.Any())
                {
                    var tasks = batch.Select(DownloadFile).ToArray();
                    Task.WaitAll(tasks);
                }
                Thread.Sleep(200);
            }
        }

        public static int downloaded = 0;
        private static async Task DownloadFile(Sequence seq)
        {
            var client = new HttpClient();

            try
            {
                HttpResponseMessage responseMessage = await client.GetAsync(MakeUrl(seq.Seq));

                var result = await responseMessage.Content.ReadAsStringAsync();


                File.WriteAllText(FileName(seq), result);
                Console.WriteLine($"{seq.Id}, {Interlocked.Increment(ref downloaded)}");
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex);
            }
        }

        private static string FileName(Sequence seq)
        {
            return Path.Combine(downloadFolder, seq.Id + ".json");
        }

        static void Main()
        {
            var lines = File.ReadAllLines(inputTestCsv);
            var items = lines.Skip(1)
                .Select(l =>
                {
                    var seqStart = l.IndexOf("\"", StringComparison.Ordinal);
                    return new Sequence {Id = l.Substring(0, l.IndexOf(",", StringComparison.Ordinal)), Seq = l.Substring(seqStart + 1, l.LastIndexOf("\"", StringComparison.Ordinal) - seqStart - 1)};
                }).ToArray();

            if (!Directory.Exists(downloadFolder))
                Directory.CreateDirectory(downloadFolder);

            DownloadFiles(items);
            Console.WriteLine("ready");
            Console.ReadKey();

        }

    }
}
