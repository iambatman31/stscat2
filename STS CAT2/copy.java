import java.net.URL;
import java.util.Scanner;

public class URLContentDisplay {
    public static void main(String[] args) {
        try {
            String urlString = "https://raw.githubusercontent.com/iambatman31/sts/main/celebrity_problem.java";

            // Create a URL object and open a Scanner on it
            try (Scanner scanner = new Scanner(new URL(urlString).openStream())) {
                // Read and print the content from the URL
                while (scanner.hasNextLine()) {
                    System.out.println(scanner.nextLine());
                }
            }
        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
        }
    }
}
