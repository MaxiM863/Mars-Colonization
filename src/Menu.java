import java.util.Scanner;

public class Menu {
    public static void ShowMainMenu()
    {
        Scanner reader = new Scanner(System.in);

        System.out.println("Bienvenue dans le projet MARS-2023:");
        System.out.println("1-Calculer une trajectoire de Lambert\n");

        String menuAnswer = reader.nextLine();

        int menuChoice = Integer.parseInt(menuAnswer);

        if(menuChoice == 1)
        {
            Trajectories T= new Trajectories();
            T.OptimalSpanTimeTrajectoryLambert(2035);
        }
    }
}
