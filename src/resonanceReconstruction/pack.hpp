template<class... Ts>
class Pack;

template<class... Ts>
inline auto pack(Ts&&... ts)
{
   return Pack<std::decay_t<Ts>...>{ std::forward<Ts>(ts)... };
}

template<class... Ts>
class Pack {
   static constexpr auto indices = std::make_index_sequence<sizeof...(Ts)>{};

   template<std::size_t... Is, class... Us>
   auto add(const std::index_sequence<Is...>, const Pack<Us...> &other) const
   {
      static_assert(sizeof...(Ts) == sizeof...(Us), "");
      return pack((std::get<Is>(data) + std::get<Is>(other.data))...);
   }

   template<std::size_t... Is, class... Us>
   auto sub(const std::index_sequence<Is...>, const Pack<Us...> &other) const
   {
      static_assert(sizeof...(Ts) == sizeof...(Us), "");
      return pack((std::get<Is>(data) - std::get<Is>(other.data))...);
   }

   template<std::size_t... Is, class ScalarValue>
   auto mul(const std::index_sequence<Is...>, const ScalarValue &value) const
   {
      return pack((std::get<Is>(data) * value)...);
   }

public:

   std::tuple<Ts...> data;
   Pack(Ts... data) : data(std::make_tuple(data...)) { }

   template<class... Us>
   auto operator+(const Pack<Us...> &other) const { return add(indices,other); }
   template<class... Us>
   auto operator-(const Pack<Us...> &other) const { return sub(indices,other); }
   template<class ScalarValue>
   auto operator*(const ScalarValue &value) const { return mul(indices,value); }
};
