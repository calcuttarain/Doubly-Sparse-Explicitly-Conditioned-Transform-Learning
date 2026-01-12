function alpha = get_alpha(T_curr, D_curr, T_ant, D_ant, Y, alpha_choice)
    switch alpha_choice

        case "exact_line_search"
            alpha_eln = 0.25 * norm(D_curr, 'fro')^2 / norm(D_curr * Y, 'fro')^2;

            alpha = alpha_eln;

        case "barzilai_borwein"
            delta_D = D_curr - D_ant;
            delta_T = T_curr - T_ant;
            alpha_short_bb = trace(delta_T' * delta_D) / trace(delta_D' * delta_D);
            alpha_long_bb = trace(delta_T' * delta_T) / trace(delta_T' * delta_D);

            alpha = max([alpha_short_bb alpha_long_bb]);

        case "long_barzilai_borwein"
            delta_D = D_curr - D_ant;
            delta_T = T_curr - T_ant;
            alpha_long_bb = trace(delta_T' * delta_T) / trace(delta_T' * delta_D);

            alpha = alpha_long_bb;

        case "short_barzilai_borwein"
            delta_D = D_curr - D_ant;
            delta_T = T_curr - T_ant;
            alpha_short_bb = trace(delta_T' * delta_D) / trace(delta_D' * delta_D);

            alpha = alpha_short_bb;


        case "lipschitz_constant"
            lipschitz_constant = 2 * norm(Y, 2)^2;
            
            alpha = 1 / lipschitz_constant;
    end
end

