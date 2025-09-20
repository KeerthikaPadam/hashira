import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PolynomialSecretFinder {

    public static BigInteger baseToDecimal(String value, int base) {
        BigInteger result = BigInteger.ZERO;
        BigInteger baseBI = BigInteger.valueOf(base);
        for (int i = 0; i < value.length(); i++) {
            char digit = value.charAt(i);
            int digitValue;
            if (digit >= '0' && digit <= '9') {
                digitValue = digit - '0';
            } else {
                digitValue = Character.toLowerCase(digit) - 'a' + 10;
            }
            result = result.multiply(baseBI).add(BigInteger.valueOf(digitValue));
        }
        return result;
    }

    public static BigInteger lagrangeInterpolation(List<Point> points, int k) {
        BigInteger secret = BigInteger.ZERO;
        for (int i = 0; i < k && i < points.size(); i++) {
            BigInteger xi = points.get(i).x;
            BigInteger yi = points.get(i).y;
            BigInteger numerator = BigInteger.ONE;
            BigInteger denominator = BigInteger.ONE;
            for (int j = 0; j < k && j < points.size(); j++) {
                if (i != j) {
                    BigInteger xj = points.get(j).x;
                    numerator = numerator.multiply(BigInteger.ZERO.subtract(xj));
                    denominator = denominator.multiply(xi.subtract(xj));
                }
            }
            secret = secret.add(yi.multiply(numerator).divide(denominator));
        }
        return secret;
    }

    static class Point {
        BigInteger x, y;
        Point(BigInteger x, BigInteger y) {
            this.x = x;
            this.y = y;
        }
    }

    public static TestCase parseJSON(String jsonContent) {
        TestCase testCase = new TestCase();
        Pattern nPattern = Pattern.compile("\"n\"\\s*:\\s*(\\d+)");
        Pattern kPattern = Pattern.compile("\"k\"\\s*:\\s*(\\d+)");
        Matcher nMatcher = nPattern.matcher(jsonContent);
        Matcher kMatcher = kPattern.matcher(jsonContent);
        if (nMatcher.find()) testCase.n = Integer.parseInt(nMatcher.group(1));
        if (kMatcher.find()) testCase.k = Integer.parseInt(kMatcher.group(1));
        Pattern rootPattern = Pattern.compile("\"(\\d+)\"\\s*:\\s*\\{[^}]*\"base\"\\s*:\\s*\"(\\d+)\"[^}]*\"value\"\\s*:\\s*\"([^\"]+)\"[^}]*\\}");
        Matcher rootMatcher = rootPattern.matcher(jsonContent);
        while (rootMatcher.find()) {
            int rootId = Integer.parseInt(rootMatcher.group(1));
            int base = Integer.parseInt(rootMatcher.group(2));
            String value = rootMatcher.group(3);
            RootEntry root = new RootEntry();
            root.id = rootId;
            root.base = base;
            root.value = value;
            testCase.roots.add(root);
        }
        return testCase;
    }

    static class TestCase {
        int n, k;
        List<RootEntry> roots = new ArrayList<>();
    }

    static class RootEntry {
        int id, base;
        String value;
    }

    public static String readFile(String filename) throws IOException {
        StringBuilder content = new StringBuilder();
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = reader.readLine()) != null) {
                content.append(line).append("\n");
            }
        }
        return content.toString();
    }

    public static void printSecret(String filename) {
        try {
            String jsonContent = readFile(filename);
            TestCase testCase = parseJSON(jsonContent);
            List<Point> points = new ArrayList<>();
            for (RootEntry root : testCase.roots) {
                BigInteger x = BigInteger.valueOf(root.id);
                BigInteger y = baseToDecimal(root.value, root.base);
                points.add(new Point(x, y));
            }
            BigInteger secret = lagrangeInterpolation(points, testCase.k);
            System.out.println(secret);  // Only prints the constant
        } catch (IOException e) {
            // Print nothing if error
        }
    }

    public static void main(String[] args) {
        if (args.length >= 2) {
            printSecret(args[0]);
            printSecret(args[1]);
        } else {
            printSecret("input1.json");
            printSecret("input2.json");
        }
    }
}
