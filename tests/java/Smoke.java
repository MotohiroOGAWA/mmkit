import org.mmkit.*;

public class Smoke {
    public static void main(String[] args) {
        System.loadLibrary("mmkit_java");
        try (Compound compound = Compound.from_smiles("CCO");
             Formula formula = compound.formula();
             Adduct adduct = Adduct.parse("[M+H]+")) {
            if (!formula.to_string().equals("C2H6O")) throw new AssertionError("formula");
            if (Math.abs(adduct.apply_to_mz(compound.exact_mass()) - 47.049689812) > 1e-6)
                throw new AssertionError("m/z");
            try {
                Adduct.parse("bad");
                throw new AssertionError("invalid input accepted");
            } catch (IllegalArgumentException expected) { }
            System.out.println("Java JNI OK: " + formula.to_string());
        }
    }
}
