// Node.js (ES2020+) — Read JSON from stdin, reconstruct f(0) via Lagrange interpolation with BigInt

// BigInt GCD (Euclidean algorithm)
function bigGcd(a, b) {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    const t = b;
    b = a % b;
    a = t;
  }
  return a;
}

// Rational number with BigInt numerator/denominator
class Fraction {
  constructor(n, d = 1n) {
    if (d === 0n) throw new Error("Zero denominator");
    // normalize sign to denominator > 0
    if (d < 0n) {
      n = -n;
      d = -d;
    }
    const g = bigGcd(n, d);
    this.n = n / g;
    this.d = d / g;
  }
  static fromBigInt(x) {
    return new Fraction(x, 1n);
  }
  add(other) {
    const n = this.n * other.d + other.n * this.d;
    const d = this.d * other.d;
    return new Fraction(n, d);
  }
  mul(other) {
    const n = this.n * other.n;
    const d = this.d * other.d;
    return new Fraction(n, d);
  }
  toBigIntExact() {
    if (this.n % this.d !== 0n) {
      throw new Error(`Non-integer result: ${this.n}/${this.d}`);
    }
    return this.n / this.d;
  }
}

// Parse a string in base 2..36 to BigInt
function parseInBase(str, base) {
  const b = BigInt(base);
  let res = 0n;
  for (const ch of str.trim().toLowerCase()) {
    let digit;
    if (ch >= '0' && ch <= '9') digit = BigInt(ch.charCodeAt(0) - 48);
    else if (ch >= 'a' && ch <= 'z') digit = BigInt(10 + ch.charCodeAt(0) - 97);
    else if (ch === ' ') continue;
    else throw new Error(`Invalid digit: ${ch}`);
    if (digit >= b) throw new Error(`Digit ${ch} out of range for base ${base}`);
    res = res * b + digit;
  }
  return res;
}

// Build points [(x_i, y_i)] from JSON
function pointsFromJson(obj) {
  const { n, k } = obj.keys;
  const entries = Object.entries(obj)
    .filter(([key]) => key !== "keys")
    .map(([key, val]) => ({ x: BigInt(key), y: parseInBase(val.value, Number(val.base)) }));
  // sort by x to make subset selection deterministic
  entries.sort((a, b) => (a.x < b.x ? -1 : a.x > b.x ? 1 : 0));
  // proceed even if entries.length !== n (some judges ignore n)
  return { k, points: entries };
}

// Lagrange interpolation at x = 0 using exactly k points (first k by sorted x)
function reconstructAtZero(points, k) {
  const P = points.slice(0, k);
  let sum = Fraction.fromBigInt(0n);
  for (let i = 0; i < P.length; i++) {
    const xi = P[i].x;
    const yi = P[i].y;
    // Compute basis weight w_i(0) = Π_{j≠i} (-x_j)/(x_i - x_j)
    let num = 1n;
    let den = 1n;
    for (let j = 0; j < P.length; j++) {
      if (j === i) continue;
      const xj = P[j].x;
      num *= -xj;
      den *= (xi - xj);
      // Reduce intermittently to keep numbers smaller
      const g = bigGcd(num, den);
      if (g > 1n) {
        num /= g;
        den /= g;
      }
    }
    const term = new Fraction(yi * num, den);
    sum = sum.add(term);
  }
  return sum.toBigIntExact();
}

// Main: read stdin, parse JSON, output f(0) in the specified format
function main() {
  const fs = require('fs');
  const input = fs.readFileSync(0, 'utf8').trim();
  if (!input) return;
  
  const obj = JSON.parse(input);
  const { k, points } = pointsFromJson(obj);

  // Ensure there are enough points to proceed
  if (points.length < k) {
    console.error(`Error: Interpolation requires ${k} points, but only ${points.length} were provided.`);
    return;
  }

  const polynomialDegree = k - 1;
  const secret = reconstructAtZero(points, k);

  // --- Formatted Output ---
  console.log(`--- Processing: Test Case ---`);
  console.log(`Polynomial degree m = k-1 = ${polynomialDegree}`);
  console.log(`Using ${k} points for interpolation.`);
  console.log(``); // For the blank line
  console.log(`--- RESULT ---`);
  console.log(`The constant term 'c' is: ${secret.toString()}`);
  console.log(`--------------`);
}


if (require.main === module) {
  main();
}
